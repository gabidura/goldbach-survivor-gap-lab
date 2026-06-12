#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = uint64_t;
using u128 = __uint128_t;

static std::string now_iso() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

static double seconds_since(const std::chrono::steady_clock::time_point& t0) {
    using namespace std::chrono;
    return duration_cast<duration<double>>(steady_clock::now() - t0).count();
}

static u64 parse_u64(const std::string& s) {
    if (s.empty()) throw std::runtime_error("empty integer argument");
    u64 value = 0;
    for (char c : s) {
        if (c < '0' || c > '9') throw std::runtime_error("invalid integer: " + s);
        unsigned digit = static_cast<unsigned>(c - '0');
        if (value > (std::numeric_limits<u64>::max() - digit) / 10) {
            throw std::runtime_error("integer too large for uint64_t: " + s);
        }
        value = value * 10 + digit;
    }
    return value;
}

static u64 floor_sqrt_u64(u64 n) {
    long double r = std::sqrt(static_cast<long double>(n));
    u64 x = static_cast<u64>(r);
    while ((u128)(x + 1) * (u128)(x + 1) <= n) ++x;
    while ((u128)x * (u128)x > n) --x;
    return x;
}

class OddCompositeBitset {
public:
    OddCompositeBitset() = default;
    explicit OddCompositeBitset(u64 limit) { reset(limit); }

    void reset(u64 limit_) {
        limit = limit_;
        // One bit per odd number n, indexed by n >> 1.
        u64 bits = (limit >> 1) + 1;
        words.assign((bits + 63) >> 6, 0ULL);
    }

    void set_composite_odd(u64 n) {
        u64 idx = n >> 1;
        words[idx >> 6] |= (1ULL << (idx & 63));
    }

    bool is_composite_odd(u64 n) const {
        u64 idx = n >> 1;
        return (words[idx >> 6] >> (idx & 63)) & 1ULL;
    }

    u64 max_limit() const { return limit; }
    size_t bytes() const { return words.size() * sizeof(u64); }

private:
    u64 limit = 0;
    std::vector<u64> words;
};

class DensePrimeTable {
public:
    DensePrimeTable() = default;

    void build(u64 limit_) {
        limit = limit_;
        comp.reset(limit);
        if (limit < 3) return;
        for (u64 p = 3; p <= limit / p; p += 2) {
            if (!comp.is_composite_odd(p)) {
                u64 step = 2 * p;
                for (u64 m = p * p; m <= limit; m += step) {
                    comp.set_composite_odd(m);
                }
            }
        }
    }

    bool contains(u64 n) const { return n <= limit; }

    bool is_prime(u64 n) const {
        if (n < 2) return false;
        if (n == 2) return true;
        if ((n & 1ULL) == 0) return false;
        if (n > limit) throw std::runtime_error("DensePrimeTable queried above limit");
        return !comp.is_composite_odd(n);
    }

    u64 max_limit() const { return limit; }
    size_t bytes() const { return comp.bytes(); }

private:
    u64 limit = 0;
    OddCompositeBitset comp;
};

static u64 mod_pow_u64(u64 a, u64 d, u64 mod) {
    u64 result = 1;
    while (d) {
        if (d & 1ULL) result = static_cast<u64>((u128)result * a % mod);
        a = static_cast<u64>((u128)a * a % mod);
        d >>= 1ULL;
    }
    return result;
}

static bool is_prime_miller_u64(u64 n) {
    if (n < 2) return false;
    static const u64 small_primes[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (u64 p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    u64 d = n - 1;
    unsigned s = 0;
    while ((d & 1ULL) == 0) { d >>= 1ULL; ++s; }

    // Deterministic bases for testing all uint64_t integers.
    static const u64 bases[] = {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};
    for (u64 a : bases) {
        if (a % n == 0) continue;
        u64 x = mod_pow_u64(a % n, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (unsigned r = 1; r < s; ++r) {
            x = static_cast<u64>((u128)x * x % n);
            if (x == n - 1) { witness = false; break; }
        }
        if (witness) return false;
    }
    return true;
}

struct Config {
    std::string mode = "first-interior"; // first-interior, gap-sample, exact-gap, single
    std::string backend = "dense";       // dense, miller
    u64 startN = 4;
    u64 stopN = 1000000;
    u64 step = 2;
    std::string output_prefix = "goldbach_v03";
    bool write_full = false;
    bool check_boundary = false;
    u64 progress_every = 10000000;
};

static void print_usage() {
    std::cout <<
R"USAGE(Goldbach Lab v03

Usage examples:
  goldbach_lab_v03 --mode first-interior --maxN 10000000000 --backend dense --output-prefix goldbach_to_10B_v03
  goldbach_lab_v03 --mode gap-sample --maxN 10000000000 --step 1000000 --backend dense --output-prefix goldbach_gap_sample_to_10B_v03
  goldbach_lab_v03 --mode exact-gap --maxN 1000000 --backend dense --output-prefix goldbach_exact_gap_to_1M_v03
  goldbach_lab_v03 --mode single --N 1000000000000 --backend miller --output-prefix single_1T_v03

Options:
  --mode <first-interior|gap-sample|exact-gap|single>
  --backend <dense|miller>
  --maxN <integer>             Sets stopN. startN defaults to 4.
  --startN <integer>
  --stopN <integer>
  --N <integer>                For --mode single.
  --step <integer>             For sampled modes. Will be rounded to an even step when scanning even N.
  --output-prefix <prefix>
  --write-full <0|1>           Write every even N row. Not recommended for huge scans.
  --check-boundary <0|1>       Compute boundary-channel positivity. Slower, but useful for reports.
  --progress-every <integer>   Print progress every this many even N values. 0 disables progress.

Notes:
  dense backend stores an odd-only prime table up to stopN, using approximately stopN/16 bytes.
  miller backend uses deterministic 64-bit Miller--Rabin and stores no global prime table.
  All arithmetic is uint64_t; use values <= 18446744073709551615.
)USAGE";
}

static Config parse_args(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string key = argv[i];
        auto need = [&](const std::string& k)->std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value after " + k);
            return argv[++i];
        };
        if (key == "--help" || key == "-h") { print_usage(); std::exit(0); }
        else if (key == "--mode") cfg.mode = need(key);
        else if (key == "--backend") cfg.backend = need(key);
        else if (key == "--maxN") cfg.stopN = parse_u64(need(key));
        else if (key == "--startN") cfg.startN = parse_u64(need(key));
        else if (key == "--stopN") cfg.stopN = parse_u64(need(key));
        else if (key == "--N") { cfg.startN = cfg.stopN = parse_u64(need(key)); cfg.mode = "single"; }
        else if (key == "--step") cfg.step = parse_u64(need(key));
        else if (key == "--output-prefix") cfg.output_prefix = need(key);
        else if (key == "--write-full") cfg.write_full = (parse_u64(need(key)) != 0);
        else if (key == "--check-boundary") cfg.check_boundary = (parse_u64(need(key)) != 0);
        else if (key == "--progress-every") cfg.progress_every = parse_u64(need(key));
        else throw std::runtime_error("unknown option: " + key);
    }
    if (cfg.mode != "first-interior" && cfg.mode != "gap-sample" && cfg.mode != "exact-gap" && cfg.mode != "single")
        throw std::runtime_error("unknown mode: " + cfg.mode);
    if (cfg.backend != "dense" && cfg.backend != "miller")
        throw std::runtime_error("unknown backend: " + cfg.backend);
    if (cfg.startN > cfg.stopN) throw std::runtime_error("startN exceeds stopN");
    if (cfg.step == 0) cfg.step = 2;
    if (cfg.step & 1ULL) ++cfg.step;
    if (cfg.startN < 4) cfg.startN = 4;
    if (cfg.startN & 1ULL) ++cfg.startN;
    if (cfg.stopN & 1ULL) --cfg.stopN;
    return cfg;
}

class PrimeOracle {
public:
    explicit PrimeOracle(std::string backend_) : backend(std::move(backend_)) {}

    void build(u64 limit) {
        if (backend == "dense") dense.build(limit);
    }

    bool is_prime(u64 n) const {
        if (backend == "dense") {
            if (dense.contains(n)) return dense.is_prime(n);
            return is_prime_miller_u64(n);
        }
        return is_prime_miller_u64(n);
    }

    size_t bytes() const { return backend == "dense" ? dense.bytes() : 0; }
    u64 dense_limit() const { return backend == "dense" ? dense.max_limit() : 0; }

private:
    std::string backend;
    DensePrimeTable dense;
};

struct PairResult {
    bool found = false;
    u64 N = 0;
    u64 z = 0;
    u64 x = 0;
    u64 y = 0;
    u64 offset = 0;
    long double ratio = 0.0L;
};

static PairResult first_interior_pair(u64 N, const PrimeOracle& P) {
    PairResult r;
    r.N = N;
    r.z = floor_sqrt_u64(N);
    if (N < 6) return r;
    u64 start = r.z + 1;
    if (start <= 2) start = 3;
    if ((start & 1ULL) == 0) ++start;
    u64 end = N - r.z - 1; // x < N-z
    for (u64 x = start; x <= end; x += 2) {
        if (P.is_prime(x)) {
            u64 y = N - x;
            if (y > r.z && P.is_prime(y)) {
                r.found = true;
                r.x = x;
                r.y = y;
                r.offset = x - r.z;
                r.ratio = static_cast<long double>(x) / static_cast<long double>(r.z == 0 ? 1 : r.z);
                return r;
            }
        }
        if (x > std::numeric_limits<u64>::max() - 2) break;
    }
    return r;
}

static bool boundary_positive(u64 N, u64 z, const PrimeOracle& P) {
    if (z >= 2 && N >= 4 && P.is_prime(N - 2)) return true;
    for (u64 q = 3; q <= z; q += 2) {
        if (P.is_prime(q) && q < N && P.is_prime(N - q)) return true;
    }
    return false;
}

struct GapResult {
    u64 N = 0;
    u64 z = 0;
    u64 max_empty_run = 0;
    long double ratio_to_N = 0.0L;
    u64 interior_count = 0;
    u64 first_x = 0;
    u64 last_x = 0;
};

static GapResult complete_window_gap(u64 N, const PrimeOracle& P) {
    GapResult g;
    g.N = N;
    g.z = floor_sqrt_u64(N);
    if (N <= 4 || N <= 2 * g.z + 1) {
        g.max_empty_run = 0;
        return g;
    }
    u64 start = g.z + 1;
    u64 end = N - g.z - 1;
    u64 prev = 0;
    bool have_prev = false;
    g.max_empty_run = 0;

    // Prime x must be odd for even N>4. Still compute empty runs over all integer positions.
    u64 x0 = start;
    if (x0 <= 2) x0 = 3;
    if ((x0 & 1ULL) == 0) ++x0;

    for (u64 x = x0; x <= end; x += 2) {
        if (P.is_prime(x) && P.is_prime(N - x)) {
            if (!have_prev) {
                if (x > start) g.max_empty_run = std::max(g.max_empty_run, x - start);
                g.first_x = x;
            } else {
                if (x > prev + 1) g.max_empty_run = std::max(g.max_empty_run, x - prev - 1);
            }
            prev = x;
            have_prev = true;
            ++g.interior_count;
            g.last_x = x;
        }
        if (x > std::numeric_limits<u64>::max() - 2) break;
    }
    if (!have_prev) {
        g.max_empty_run = (end >= start) ? (end - start + 1) : 0;
    } else {
        if (end > prev) g.max_empty_run = std::max(g.max_empty_run, end - prev);
    }
    g.ratio_to_N = static_cast<long double>(g.max_empty_run) / static_cast<long double>(N == 0 ? 1 : N);
    return g;
}

static std::string fmt_ld(long double x, int precision = 12) {
    std::ostringstream oss;
    oss << std::setprecision(precision) << static_cast<long double>(x);
    return oss.str();
}

static void run_first_interior(const Config& cfg, const PrimeOracle& P, double sieve_seconds) {
    auto t0 = std::chrono::steady_clock::now();
    const std::string prefix = cfg.output_prefix;
    std::ofstream records(prefix + "_first_interior_records.csv");
    records << "record_type,N,z,first_x,first_y,offset,ratio\n";
    std::ofstream failures(prefix + "_failures.csv");
    failures << "N,z,boundary_positive,classification\n";
    std::ofstream full;
    if (cfg.write_full) {
        full.open(prefix + "_full_first_scan.csv");
        full << "N,z,first_x,first_y,offset,ratio,boundary_positive,classification\n";
    }

    u64 tested = 0, interior_positive = 0, boundary_pos_count = 0, boundary_only = 0, no_rep = 0;
    PairResult worst_x, worst_offset, worst_ratio;
    bool have_worst_x = false, have_worst_offset = false, have_worst_ratio = false;

    for (u64 N = cfg.startN; N <= cfg.stopN; N += 2) {
        ++tested;
        PairResult r = first_interior_pair(N, P);
        bool bpos = false;
        std::string cls;
        if (r.found) {
            ++interior_positive;
            cls = "interior_positive";
            if (cfg.check_boundary) {
                bpos = boundary_positive(N, r.z, P);
                if (bpos) ++boundary_pos_count;
            }
            if (!have_worst_x || r.x > worst_x.x) {
                worst_x = r; have_worst_x = true;
                records << "first_x," << N << ',' << r.z << ',' << r.x << ',' << r.y << ',' << r.offset << ',' << fmt_ld(r.ratio) << "\n";
            }
            if (!have_worst_offset || r.offset > worst_offset.offset) {
                worst_offset = r; have_worst_offset = true;
                records << "offset," << N << ',' << r.z << ',' << r.x << ',' << r.y << ',' << r.offset << ',' << fmt_ld(r.ratio) << "\n";
            }
            if (!have_worst_ratio || r.ratio > worst_ratio.ratio) {
                worst_ratio = r; have_worst_ratio = true;
                records << "ratio," << N << ',' << r.z << ',' << r.x << ',' << r.y << ',' << r.offset << ',' << fmt_ld(r.ratio) << "\n";
            }
        } else {
            u64 z = floor_sqrt_u64(N);
            bpos = boundary_positive(N, z, P);
            if (bpos) { ++boundary_pos_count; ++boundary_only; cls = "boundary_only"; }
            else { ++no_rep; cls = "no_representation_found"; }
            failures << N << ',' << z << ',' << (bpos ? 1 : 0) << ',' << cls << "\n";
        }
        if (cfg.write_full) {
            full << N << ',' << r.z << ',' << r.x << ',' << r.y << ',' << r.offset << ',' << fmt_ld(r.ratio) << ',' << (bpos ? 1 : 0) << ',' << cls << "\n";
        }
        if (cfg.progress_every && tested % cfg.progress_every == 0) {
            std::cerr << "tested " << tested << " even N, current N=" << N << ", elapsed=" << seconds_since(t0) << "s\n";
        }
        if (N > std::numeric_limits<u64>::max() - 2) break;
    }
    double scan_seconds = seconds_since(t0);

    std::ofstream report(prefix + "_report.txt");
    report << "Goldbach Lab v03 first-interior scan\n";
    report << "====================================\n\n";
    report << "Generated: " << now_iso() << "\n";
    report << "Mode: " << cfg.mode << "\n";
    report << "Backend: " << cfg.backend << "\n";
    report << "Start N: " << cfg.startN << "\n";
    report << "Stop N: " << cfg.stopN << "\n";
    report << "Even N tested: " << tested << "\n";
    report << "Interior positive: " << interior_positive << "\n";
    if (cfg.check_boundary) report << "Boundary positive: " << boundary_pos_count << "\n";
    else report << "Boundary positive: not computed except failures\n";
    report << "Boundary-only: " << boundary_only << "\n";
    report << "No representation found: " << no_rep << "\n\n";
    if (have_worst_x) {
        report << "Worst first interior x case: N=" << worst_x.N << ", z=" << worst_x.z << ", first pair=(" << worst_x.x << ',' << worst_x.y << "), offset=" << worst_x.offset << ", ratio=" << fmt_ld(worst_x.ratio) << "\n";
    }
    if (have_worst_offset) {
        report << "Worst first-x offset above sqrt(N): N=" << worst_offset.N << ", z=" << worst_offset.z << ", first pair=(" << worst_offset.x << ',' << worst_offset.y << "), offset=" << worst_offset.offset << ", ratio=" << fmt_ld(worst_offset.ratio) << "\n";
    }
    if (have_worst_ratio) {
        report << "Worst first_x/z ratio: N=" << worst_ratio.N << ", z=" << worst_ratio.z << ", first pair=(" << worst_ratio.x << ',' << worst_ratio.y << "), offset=" << worst_ratio.offset << ", ratio=" << fmt_ld(worst_ratio.ratio) << "\n";
    }
    report << "\nOutput files:\n";
    report << prefix << "_first_interior_records.csv\n";
    report << prefix << "_failures.csv\n";
    if (cfg.write_full) report << prefix << "_full_first_scan.csv\n";
    report << "\nTiming seconds:\n";
    report << "sieve/build: " << sieve_seconds << "\n";
    report << "scan: " << scan_seconds << "\n";
    report << "total: " << (sieve_seconds + scan_seconds) << "\n";
    report << "\nNote: This scan tests first complete-sieve interior Goldbach pairs x>N^0.5 and N-x>N^0.5.\n";

    std::cout << "Wrote " << prefix << "_report.txt\n";
}

static void run_gap_scan(const Config& cfg, const PrimeOracle& P, double sieve_seconds) {
    auto t0 = std::chrono::steady_clock::now();
    const std::string prefix = cfg.output_prefix;
    bool sampled = (cfg.mode == "gap-sample");
    std::ofstream rows(prefix + (sampled ? "_gap_sample.csv" : "_exact_gap_records.csv"));
    rows << "N,z,max_empty_run,ratio_to_N,interior_count,first_x,last_x\n";

    u64 tested = 0;
    GapResult worst_gap, worst_ratio;
    bool have_gap = false, have_ratio = false;

    for (u64 N = cfg.startN; N <= cfg.stopN; N += (sampled ? cfg.step : 2)) {
        if (N & 1ULL) ++N;
        ++tested;
        GapResult g = complete_window_gap(N, P);
        bool write_row = sampled;
        if (!have_gap || g.max_empty_run > worst_gap.max_empty_run) { worst_gap = g; have_gap = true; write_row = true; }
        if (!have_ratio || g.ratio_to_N > worst_ratio.ratio_to_N) { worst_ratio = g; have_ratio = true; write_row = true; }
        if (write_row) {
            rows << g.N << ',' << g.z << ',' << g.max_empty_run << ',' << fmt_ld(g.ratio_to_N) << ',' << g.interior_count << ',' << g.first_x << ',' << g.last_x << "\n";
        }
        if (cfg.progress_every && tested % cfg.progress_every == 0) {
            std::cerr << "gap-tested " << tested << " values, current N=" << N << ", elapsed=" << seconds_since(t0) << "s\n";
        }
        if (N > std::numeric_limits<u64>::max() - cfg.step) break;
    }
    double scan_seconds = seconds_since(t0);

    std::ofstream report(prefix + "_report.txt");
    report << "Goldbach Lab v03 " << (sampled ? "sampled" : "exact") << " complete-window gap scan\n";
    report << "==================================================\n\n";
    report << "Generated: " << now_iso() << "\n";
    report << "Mode: " << cfg.mode << "\n";
    report << "Backend: " << cfg.backend << "\n";
    report << "Start N: " << cfg.startN << "\n";
    report << "Stop N: " << cfg.stopN << "\n";
    if (sampled) report << "Step: " << cfg.step << "\n";
    report << "Values tested: " << tested << "\n\n";
    if (have_gap) {
        report << "Worst max_empty_run: N=" << worst_gap.N << ", z=" << worst_gap.z << ", max_empty_run=" << worst_gap.max_empty_run << ", ratio_to_N=" << fmt_ld(worst_gap.ratio_to_N) << ", interior_count=" << worst_gap.interior_count << "\n";
    }
    if (have_ratio) {
        report << "Worst ratio max_empty_run/N: N=" << worst_ratio.N << ", z=" << worst_ratio.z << ", max_empty_run=" << worst_ratio.max_empty_run << ", ratio_to_N=" << fmt_ld(worst_ratio.ratio_to_N) << ", interior_count=" << worst_ratio.interior_count << "\n";
    }
    report << "\nOutput files:\n";
    report << prefix << (sampled ? "_gap_sample.csv" : "_exact_gap_records.csv") << "\n";
    report << "\nTiming seconds:\n";
    report << "sieve/build: " << sieve_seconds << "\n";
    report << "scan: " << scan_seconds << "\n";
    report << "total: " << (sieve_seconds + scan_seconds) << "\n";
    std::cout << "Wrote " << prefix << "_report.txt\n";
}

static void run_single(const Config& cfg, const PrimeOracle& P, double sieve_seconds) {
    auto t0 = std::chrono::steady_clock::now();
    u64 N = cfg.startN;
    if (N & 1ULL) throw std::runtime_error("single N must be even");
    PairResult r = first_interior_pair(N, P);
    bool bpos = false;
    if (!r.found || cfg.check_boundary) bpos = boundary_positive(N, floor_sqrt_u64(N), P);
    GapResult g = complete_window_gap(N, P);
    double scan_seconds = seconds_since(t0);

    std::ofstream report(cfg.output_prefix + "_report.txt");
    report << "Goldbach Lab v03 single-N report\n";
    report << "================================\n\n";
    report << "Generated: " << now_iso() << "\n";
    report << "Backend: " << cfg.backend << "\n";
    report << "N: " << N << "\n";
    report << "z=floor(sqrt(N)): " << floor_sqrt_u64(N) << "\n";
    report << "Interior positive: " << (r.found ? 1 : 0) << "\n";
    if (r.found) report << "First interior pair: (" << r.x << ',' << r.y << "), offset=" << r.offset << ", ratio=" << fmt_ld(r.ratio) << "\n";
    report << "Boundary positive: " << (bpos ? 1 : 0) << "\n";
    report << "Complete-window max_empty_run: " << g.max_empty_run << "\n";
    report << "Complete-window ratio_to_N: " << fmt_ld(g.ratio_to_N) << "\n";
    report << "Interior count: " << g.interior_count << "\n";
    report << "Timing seconds: sieve/build=" << sieve_seconds << ", scan=" << scan_seconds << ", total=" << (sieve_seconds + scan_seconds) << "\n";
    std::cout << "Wrote " << cfg.output_prefix << "_report.txt\n";
}

int main(int argc, char** argv) {
    try {
        Config cfg = parse_args(argc, argv);
        std::cerr << "Goldbach Lab v03\n";
        std::cerr << "mode=" << cfg.mode << ", backend=" << cfg.backend << ", startN=" << cfg.startN << ", stopN=" << cfg.stopN << "\n";

        auto t_sieve = std::chrono::steady_clock::now();
        PrimeOracle P(cfg.backend);
        if (cfg.backend == "dense") {
            long double mem_bytes = static_cast<long double>(cfg.stopN) / 16.0L;
            std::cerr << "estimated dense odd-only bitset memory: " << std::fixed << std::setprecision(2) << (double)(mem_bytes / (1024.0L * 1024.0L)) << " MiB\n";
            P.build(cfg.stopN);
            std::cerr << "dense table built: " << (P.bytes() / (1024.0 * 1024.0)) << " MiB\n";
        }
        double sieve_seconds = seconds_since(t_sieve);

        if (cfg.mode == "first-interior") run_first_interior(cfg, P, sieve_seconds);
        else if (cfg.mode == "gap-sample" || cfg.mode == "exact-gap") run_gap_scan(cfg, P, sieve_seconds);
        else if (cfg.mode == "single") run_single(cfg, P, sieve_seconds);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n\n";
        print_usage();
        return 1;
    }
}
