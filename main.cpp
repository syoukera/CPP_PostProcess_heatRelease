// C++17
#include "cantera/zerodim.h"
#include <bits/stdc++.h>
#include <iostream>

using namespace std;
using namespace Cantera;

using Real = double;           // Fortran 側が double の想定
static constexpr bool needs_byteswap = false; // 必要なら true に

// バイトスワップ（4/8 バイト）
template<class T>
static inline void bswap_inplace(T& v) {
    auto* p = reinterpret_cast<unsigned char*>(&v);
    std::reverse(p, p + sizeof(T));
}

// Fortran unformatted sequentls ial の 1 レコード読み込み＋末尾マーカー検証
// 期待バイト数 expected_bytes のペイロードを data に読み込む。
// 先頭マーカーが 4B か 8B か自動判定（よくある 2 通りに対応）。
static void read_fortran_record(ifstream& ifs, char* data, size_t expected_bytes) {
    if (!ifs) throw runtime_error("stream not open");

    auto read_len = [&](int64_t& L, size_t nbytes) {
        ifs.read(reinterpret_cast<char*>(&L), nbytes);
        if (!ifs) throw runtime_error("failed to read record header");
        if (needs_byteswap) bswap_inplace(L);
    };

    // まず 4B ヘッダを読んでみる
    streampos pos0 = ifs.tellg();
    int32_t L4 = 0;
    ifs.read(reinterpret_cast<char*>(&L4), 4);
    if (!ifs) throw runtime_error("failed to read record header (4B)");
    if (needs_byteswap) bswap_inplace(L4);

    bool ok = false;
    if (static_cast<uint32_t>(L4) == expected_bytes) {
        // 4B マーカー方式
        ifs.read(data, expected_bytes);
        if (!ifs) throw runtime_error("failed to read record payload (4B)");
        int32_t tail = 0;
        ifs.read(reinterpret_cast<char*>(&tail), 4);
        if (!ifs) throw runtime_error("failed to read record trailer (4B)");
        if (needs_byteswap) bswap_inplace(tail);
        if (tail != L4) throw runtime_error("record trailer mismatch (4B)");
        ok = true;
    } else {
        // 8B の可能性を試す
        ifs.clear();
        ifs.seekg(pos0);
        int64_t L8 = 0;
        read_len(L8, 8);
        if (static_cast<uint64_t>(L8) != expected_bytes) {
            throw runtime_error("record length mismatch (neither 4B nor 8B)");
        }
        ifs.read(data, expected_bytes);
        if (!ifs) throw runtime_error("failed to read record payload (8B)");
        int64_t tail = 0;
        read_len(tail, 8);
        if (tail != L8) throw runtime_error("record trailer mismatch (8B)");
        ok = true;
    }
    if (!ok) throw runtime_error("unknown record format");
}

// 0 埋め文字列
static string zero_pad(long long value, int width) {
    ostringstream oss;
    oss << setw(width) << setfill('0') << value;
    return oss.str();
}

// 区間構造体（両端含む, 1-based 入力を前提）
struct Range { int sta, end; };

// 区間交差判定
static inline bool disjoint(const Range& a, const Range& b) {
    return (a.end < b.sta) || (a.sta > b.end);
}
// static inline int clamp(int v, int lo, int hi) { return std::max(lo, std::min(v, hi)); }

// 線形化 index（Fortran と同じ i 最速）
// i in [i0, i0+ni-1], j in [j0, ...], k in [k0, ...]
static inline size_t idx3(int i, int j, int k,
                          int i0, int j0, int k0,
                          int ni, int nj, int /*nk*/) {
    return size_t(i - i0) + size_t(ni) * (size_t(j - j0) + size_t(nj) * size_t(k - k0));
}
static inline size_t idx4(int i, int j, int k, int s,
                          int i0, int j0, int k0, int s0,
                          int ni, int nj, int nk, int /*ns*/) {
    // i + ni*( j + nj*( k + nk*(s) ) ), s が最も遅い（Fortran と一致）
    return size_t(i - i0)
         + size_t(ni) * ( size_t(j - j0)
         + size_t(nj) * ( size_t(k - k0)
         + size_t(nk) * size_t(s - s0) ) );
}

// チャンネル名を定義
static const vector<string> channelNames = {
    "HeatRelease [J/m3/s]"
};

vector<Real> read_coord(const string& fname, int n, int bd) {
    vector<Real> g(n + 2*bd);
    ifstream ifs(fname);
    if (!ifs) throw runtime_error("Cannot open " + fname);
    for (int i = -bd; i < n+bd; i++) {
        int idx; Real u, gi, du, dgi;
        ifs >> idx >> u >> gi >> du >> dgi;
        g[i+bd] = gi;
    }
    return g;
}

// Base64エンコード
static const char b64_table[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

std::string base64_encode(const unsigned char* data, size_t len) {
    std::string out;
    out.reserve((len + 2) / 3 * 4);
    for (size_t i = 0; i < len; i += 3) {
        uint32_t v = data[i] << 16;
        if (i + 1 < len) v |= data[i+1] << 8;
        if (i + 2 < len) v |= data[i+2];
        out.push_back(b64_table[(v >> 18) & 0x3F]);
        out.push_back(b64_table[(v >> 12) & 0x3F]);
        if (i + 1 < len) out.push_back(b64_table[(v >> 6) & 0x3F]);
        else out.push_back('=');
        if (i + 2 < len) out.push_back(b64_table[v & 0x3F]);
        else out.push_back('=');
    }
    return out;
}

void write_vts_binary(const string& filename,
                      const vector<Real>& sf,
                      int NX, int NY, int NZ,
                      int x0, int y0, int z0,
                      const vector<Real>& xg,
                      const vector<Real>& yg,
                      const vector<Real>& zg,
                      int nch)
{
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    ofs << R"(<?xml version="1.0"?>)" << "\n";
    ofs << R"(<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">)" << "\n";

    ofs << "  <StructuredGrid WholeExtent=\""
        << "0 " << NX-1 << " "
        << "0 " << NY-1 << " "
        << "0 " << NZ-1 << "\">\n";

    ofs << "    <Piece Extent=\""
        << "0 " << NX-1 << " "
        << "0 " << NY-1 << " "
        << "0 " << NZ-1 << "\">\n";

    // --- 座標 ---
    ofs << "      <Points>\n";
    ofs << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">\n";

    {
        size_t npts = NX * NY * NZ;
        vector<double> coords(3 * npts);

        size_t idx = 0;
        for (int k = 0; k < NZ; ++k) {
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    coords[idx++] = xg[x0-1 + i];
                    coords[idx++] = yg[y0-1 + j];
                    coords[idx++] = zg[z0-1 + k];
                }
            }
        }

        // バイト数を先頭に付加
        uint32_t nbytes = static_cast<uint32_t>(coords.size() * sizeof(double));
        std::ostringstream raw;
        raw.write(reinterpret_cast<const char*>(&nbytes), sizeof(uint32_t));
        raw.write(reinterpret_cast<const char*>(coords.data()), nbytes);

        std::string encoded = base64_encode(
            reinterpret_cast<const unsigned char*>(raw.str().data()),
            raw.str().size());

        ofs << encoded << "\n";
    }

    ofs << "        </DataArray>\n";
    ofs << "      </Points>\n";

    // --- スカラー値 ---
    ofs << "      <PointData>\n";
    size_t npts = NX * NY * NZ;
    for (int c = 0; c < nch; ++c) {
        string name = (c < (int)channelNames.size()) ? channelNames[c] : ("var" + std::to_string(c+1));

        ofs << "        <DataArray type=\"Float64\" Name=\"" << name
            << "\" NumberOfComponents=\"1\" format=\"binary\">\n";

        vector<double> vals(npts);
        for (size_t n = 0; n < npts; ++n) {
            size_t idx = n + size_t(npts) * size_t(c);
            vals[n] = sf[idx];
        }

        uint32_t nbytes = static_cast<uint32_t>(vals.size() * sizeof(double));
        std::ostringstream raw;
        raw.write(reinterpret_cast<const char*>(&nbytes), sizeof(uint32_t));
        raw.write(reinterpret_cast<const char*>(vals.data()), nbytes);

        std::string encoded = base64_encode(
            reinterpret_cast<const unsigned char*>(raw.str().data()),
            raw.str().size());

        ofs << encoded << "\n";
        ofs << "        </DataArray>\n";
    }
    ofs << "      </PointData>\n";

    ofs << "    </Piece>\n";
    ofs << "  </StructuredGrid>\n";
    ofs << "</VTKFile>\n";
}

int main() {
    // ======= 設定（必要に応じて実行時引数や設定ファイル化してください） =======
    const int nx = /* 全体 x サイズ */ 1000;
    const int ny = /* 全体 y サイズ */ 500;
    const int nz = /* 全体 z サイズ */ 40;

    const int iprocs = 25, jprocs = 20, kprocs = 4; // プロセス分割
    const int ibd = 4, jbd = 4, kbd = 4;          // ハロー幅
    const int nf  = 26;                            // y の種数（>=19 を想定）

    // 出力ウィンドウ（Fortran では 1 始まり。ここも 1 始まりで指定）
    // const Range XR{1, 2000}, YR{1, 550}, ZR{1, nz};
    const Range XR{1, 800}, YR{1, 500}, ZR{1, nz};

    // 読むステップ範囲（Fortran の step0:step2:step1 に相当）
    const int step0 = 640000, step1 = 640000, step2 = 5000; // 例：単一ステップ

    // 出力配列（ウィンドウだけ確保）。C++ は 0 始まりに直す。
    const int x0 = XR.sta; const int x1 = XR.end;  // 1-based
    const int y0 = YR.sta; const int y1 = YR.end;
    const int z0 = ZR.sta; const int z1 = ZR.end;
    const int NX = x1 - x0 + 1;
    const int NY = y1 - y0 + 1;
    const int NZ = z1 - z0 + 1;

    vector<Real> xg = read_coord("../../xg.dat", nx, ibd);
    vector<Real> yg = read_coord("../../yg.dat", ny, jbd);
    vector<Real> zg = read_coord("../../zg.dat", nz, kbd);

    // vector<Real> p_global(size_t(NX) * NY * NZ, Real(0));
    // vector<Real> sf(size_t(NX) * NY * NZ * nf, Real(0));
    
    // define cantera object
    auto sol = newSolution("Tamaoki_reduced_combine_IDT_Su.yaml", "", "none");
    auto gas = sol->thermo();
    auto kin = sol->kinetics();

    // buffefr for calculation
    // size_t nf = gas->nSpecies();
    vector<double> wdot(nf);
    vector<double> h_bar(nf);
    
    // define output vector
    vector<Real> heat_release(size_t(NX) * NY * NZ, Real(0));


    auto sfind = [&](int i, int j, int k, int c) -> size_t {
        // i,j,k は 1-based（グローバル座標）、c は 1-based チャネル番号
        // C++ 配列は 0-based に変換
        int ii = i - x0; int jj = j - y0; int kk = k - z0; int cc = c - 1;
        // i 最速, 次に j, 次に k, 最後に channel
        return size_t(ii) + size_t(NX) * ( size_t(jj) + size_t(NY) * ( size_t(kk) + size_t(NZ) * size_t(cc) ) );
    };

    auto pind = [&](int i, int j, int k) -> size_t {
        // i,j,k は 1-based（グローバル座標）
        // C++ 配列は 0-based に変換
        int ii = i - x0; int jj = j - y0; int kk = k - z0;
        // i 最速, 次に j, 次に k, 最後に channel
        return size_t(ii) + size_t(NX) * ( size_t(jj) + size_t(NY) * size_t(kk) );
    };

    // ======= ステップループ =======
    for (int step = step0; step <= step1; step += step2) {
        string filenumber = zero_pad(step, 8);

        const int nprocs = iprocs * jprocs * kprocs;
        for (int myrank = 0; myrank < nprocs; ++myrank) {
            // rank -> (i,j,k) ブロック座標（Fortran の m カウントに対応）
            int myrank_i =  myrank % iprocs;
            int myrank_j = (myrank / iprocs) % jprocs;
            int myrank_k =  myrank / (iprocs * jprocs);

            // サブドメイン範囲（1-based, 両端含む）
            int ista = nx / iprocs * myrank_i + 1;
            int iend = nx / iprocs * (myrank_i + 1);
            int jsta = ny / jprocs * myrank_j + 1;
            int jend = ny / jprocs * (myrank_j + 1);
            int ksta = nz / kprocs * myrank_k + 1;
            int kend = nz / kprocs * (myrank_k + 1);

            // ローカル（ハロー含む）配列の 1-based 範囲
            int li0 = ista - ibd, li1 = iend + ibd;
            int lj0 = jsta - jbd, lj1 = jend + jbd;
            int lk0 = ksta - kbd, lk1 = kend + kbd;

            // Range IR{li0, li1}, JR{lj0, lj1}, KR{lk0, lk1};

            // 出力ウィンドウと交差しないならスキップ
            if ( disjoint(Range{ista, iend}, XR) ||
                 disjoint(Range{jsta, jend}, YR) ||
                 disjoint(Range{ksta, kend}, ZR) ) {
                continue;
            }

            // ローカルサイズ
            int lni = li1 - li0 + 1;
            int lnj = lj1 - lj0 + 1;
            int lnk = lk1 - lk0 + 1;

            // ファイルを開く
            string cpunumber = zero_pad(myrank, 5);
            string path;
            if (nprocs != 1)
                path = string("../") + cpunumber + "/f" + filenumber + ".dat";
            else
                path = string("../f") + filenumber + ".dat";

            ifstream ifs(path, ios::binary);
            if (!ifs) {
                cerr << "Cannot open file: " << path << "\n";
                return 1;
            }
            cerr << "Reading " << path << "\n";

            // 読み込みバッファを確保（Fortran 配列の線形順序と一致させる）
            const size_t n3 = size_t(lni) * lnj * lnk;
            const size_t n4 = n3 * size_t(nf);
            vector<Real> u(n3), v(n3), w(n3), r(n3), p(n3), t(n3), h(n3);
            vector<Real> y(n4);

            auto read_vec = [&](vector<Real>& a) {
                const size_t bytes = a.size() * sizeof(Real);
                read_fortran_record(ifs, reinterpret_cast<char*>(a.data()), bytes);
                if (needs_byteswap) {
                    for (auto& val : a) bswap_inplace(val);
                }
            };

            read_vec(u); read_vec(v); read_vec(w); read_vec(r);
            read_vec(p); read_vec(t); read_vec(h);
            // y は 4 次元 (i,j,k,s) を一括で 1 レコードとして書いている想定
            read_vec(y);
            ifs.close();

            // マージ（ウィンドウ＆ハローでクリップ）
            const int ii_sta = max(li0, XR.sta);
            const int ii_end = min(li1, XR.end);
            const int jj_sta = max(lj0, YR.sta);
            const int jj_end = min(lj1, YR.end);
            const int kk_sta = max(lk0, ZR.sta);
            const int kk_end = min(lk1, ZR.end);

            for (int k = kk_sta; k <= kk_end; ++k)
            for (int j = jj_sta; j <= jj_end; ++j)
            for (int i = ii_sta; i <= ii_end; ++i) {
                size_t L3 = idx3(i, j, k, li0, lj0, lk0, lni, lnj, lnk);

                // 基本物理量
                // sf[sfind(i,j,k, 1)] = u[L3];
                // sf[sfind(i,j,k, 2)] = v[L3];
                // sf[sfind(i,j,k, 3)] = w[L3];
                // sf[sfind(i,j,k, 4)] = r[L3];
                // sf[sfind(i,j,k, 5)] = p[L3];
                // sf[sfind(i,j,k, 6)] = t[L3];
                // sf[sfind(i,j,k, 7)] = h[L3];

                // 種（Fortran では 1 始まりの添字）
                auto LY = [&](int s1based)->Real {
                    size_t L4 = idx4(i, j, k, s1based, li0, lj0, lk0, 1, lni, lnj, lnk, nf);
                    return y[L4];
                };

                // get y_local
                vector<double> y_local(nf);
                for (int s = 1; s <= nf; ++s) {
                    y_local[s-1] = LY(s);
                }

                // create a reservoir for the fuel inlet, and set to pure methane.
                gas->setState_TPY(t[L3], p[L3], y_local.data());

                kin->getNetProductionRates(wdot.data()); // kmol/m^3/s
                gas->getPartialMolarEnthalpies(h_bar.data()); // J/kmol
                
                // calcualte heat relase [J/m3/s]
                double hrr_vol = 0.0;
                for (size_t k = 0; k < nf; k++) {
                    hrr_vol -= wdot[k] * h_bar[k];
                }
                
                // asign heat_release
                heat_release[pind(i,j,k)] = hrr_vol;

            }
        }
        
        string outname = "output_step_" + filenumber + ".vts";
        write_vts_binary(outname, heat_release, NX, NY, NZ, x0, y0, z0, xg, yg, zg, 1);
        cerr << "Wrote " << outname << "\n";
    }
    
    return 0;
}