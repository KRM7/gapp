/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_DISTRIBUTION_HPP
#define GAPP_UTILITY_DISTRIBUTION_HPP

#include "bit.hpp"
#include "math.hpp"
#include <random>
#include <array>
#include <bit>
#include <optional>
#include <concepts>
#include <type_traits>
#include <limits>
#include <cmath>
#include <cstdint>

namespace gapp::rng
{
    template<typename G>
    concept prng64 = requires(G gen)
    {
        requires std::uniform_random_bit_generator<G>;

        typename G::result_type;
        requires sizeof(typename G::result_type) == sizeof(uint64_t);

        { gen() } -> std::same_as<typename G::result_type>;

        requires G::min() == std::numeric_limits<typename G::result_type>::min();
        requires G::max() == std::numeric_limits<typename G::result_type>::max();
    };


    /** @returns A random float from a uniform distribution in [0.0, 1.0). */
    template<std::same_as<float> RealType, rng::prng64 G>
    constexpr RealType generate_canonical(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
    {
        return (uint64_t{ gen() } >> (32 + detail::exponent_bits<RealType>)) *
               (std::numeric_limits<RealType>::epsilon() / 2);
    }

    /** @returns A random double from a uniform distribution in [0.0, 1.0). */
    template<std::same_as<double> RealType, rng::prng64 G>
    constexpr RealType generate_canonical(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
    {
        return (uint64_t{ gen() } >> detail::exponent_bits<RealType>) *
               (std::numeric_limits<RealType>::epsilon() / 2);
    }


    /** Generates random booleans from a uniform distribution. */
    class uniform_bool_distribution
    {
    public:
        using result_type = bool;

        constexpr uniform_bool_distribution() noexcept = default;

        template<std::uniform_random_bit_generator G>
        constexpr result_type operator()(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
        {
            if (bit_pool_ == 1)
            {
                bit_pool_ = gen() | detail::msb_mask<typename G::result_type>;
            }

            const bool bit = detail::lsb(bit_pool_);
            bit_pool_ >>= 1;

            return bit;
        }

        constexpr void reset() noexcept { bit_pool_ = 1; }

        static constexpr result_type min() noexcept { return false; }
        static constexpr result_type max() noexcept { return true; }

        friend constexpr bool operator==(const uniform_bool_distribution&, const uniform_bool_distribution&) = default;

    private:
        uint64_t bit_pool_ = 1;
    };


    /** Generates random integers from a uniform distribution in [low, high]. */
    template<std::integral T>
    class uniform_int_distribution
    {
    public:
        static_assert(sizeof(T) <= sizeof(uint64_t));

        using result_type = T;

        constexpr uniform_int_distribution(T low, T high) noexcept :
            min_(low)
        {
            GAPP_ASSERT(low <= high);

            range_ = uint64_t(high) - uint64_t(low) + 1ull;

            if (range_ == 0) return;

            partitions_ = uint64_t(-1) / range_;
            threshold_ = partitions_ * range_;
        }
        
        template<rng::prng64 G>
        constexpr result_type operator()(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
        {
            if (range_ == 0) return result_type(gen());

            for (;;)
            {
                const uint64_t value = gen();
                if (value >= threshold_)
                    continue;

                return result_type(min_ + value / partitions_);
            }

            GAPP_UNREACHABLE();
        }

        constexpr void reset() noexcept {}

        constexpr result_type min() const noexcept { return result_type(min_); }
        constexpr result_type max() const noexcept { return result_type(range_ + min_ - 1ull); }

        friend constexpr bool operator==(const uniform_int_distribution& lhs, const uniform_int_distribution& rhs) noexcept
        {
            return lhs.min_ == rhs.min_ && lhs.range_ == rhs.range_;
        }

    private:
        T min_;
        uint64_t range_;
        uint64_t partitions_ = 0;
        uint64_t threshold_  = 0;
    };


    /** Generates random floating point numbers from a uniform distribution in [low, high). */
    template<std::floating_point T>
    class uniform_real_distribution
    {
    public:
        static_assert(sizeof(T) <= sizeof(uint64_t));

        using result_type = T;

        constexpr uniform_real_distribution() noexcept = default;

        constexpr uniform_real_distribution(T low, T high) noexcept :
            min_(low), range_(high - low)
        {
            GAPP_ASSERT(low <= high);
            GAPP_ASSERT(range_ <= std::numeric_limits<T>::max());
        }

        template<rng::prng64 G>
        constexpr result_type operator()(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
        {
            return min_ + range_ * rng::generate_canonical<T>(gen);
        }

        constexpr void reset() noexcept {}

        constexpr result_type min() const noexcept { return min_; }
        constexpr result_type max() const noexcept { return range_ + min_; }

        friend constexpr bool operator==(const uniform_real_distribution&, const uniform_real_distribution&) = default;

    private:
        T min_ = 0.0;
        T range_ = 1.0;
    };


    /** Generates random floating point numbers from an exponential distribution. */
    template<std::floating_point T>
    class exponential_distribution
    {
    public:
        static_assert(sizeof(T) <= sizeof(uint64_t));

        using result_type = T;

        constexpr exponential_distribution() noexcept = default;

        explicit constexpr exponential_distribution(T lambda) noexcept :
            inv_neg_lambda_(-1 / lambda)
        {
            GAPP_ASSERT(lambda > 0.0);
        }

        template<rng::prng64 G>
        constexpr result_type operator()(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
        {
            return inv_neg_lambda_ * std::log(1 - rng::generate_canonical<T>(gen));
        }

        constexpr void reset() noexcept {}

        static constexpr result_type min() noexcept { return 0.0; }
        static constexpr result_type max() noexcept { return std::numeric_limits<T>::infinity(); }

        friend constexpr bool operator==(const exponential_distribution&, const exponential_distribution&) = default;

    private:
        T inv_neg_lambda_ = -1.0;
    };


    /** Generates random floating point number from a normal distribution. */
    template<std::floating_point T>
    class normal_distribution
    {
    public:
        static_assert(sizeof(T) <= sizeof(uint64_t));

        using result_type = T;

        constexpr normal_distribution() noexcept = default;

        constexpr normal_distribution(T mean, T stddev) noexcept :
            mean_(mean), stddev_(stddev)
        {
            GAPP_ASSERT(stddev_ > 0.0);
        }

        template<rng::prng64 G>
        constexpr result_type operator()(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
        {
            // Ziggurat method, based on:
            //  Doornik, Jurgen A. "An improved ziggurat method to generate normal random samples."
            //  University of Oxford (2005): 77

            for (;;)
            {
                const uint64_t bits = gen();

                const uint64_t sign_bit = detail::extract_bits<63, 64>(bits);
                const uint64_t i        = detail::extract_bits<56, 63>(bits);
                const uint64_t u0_bits  = detail::extract_bits<56 - detail::implicit_mantissa_bits<T>, 56>(bits);

                const T z = x[i] * u0_bits * (std::numeric_limits<T>::epsilon() / 2);

                if (z < x[i + 1])
                {
                    return stddev_ * detail::set_sign_bit(z, sign_bit) + mean_;
                }

                if (i == 0)
                {
                    constexpr T tail = x[1];
                    rng::exponential_distribution<T> exp_x(tail);
                    rng::exponential_distribution<T> exp_y(1.0);

                    for (;;)
                    {
                        const T ex = exp_x(gen);
                        const T ey = exp_y(gen);
                        if (ey + ey <= ex * ex)
                            continue;

                        return stddev_ * detail::set_sign_bit(ex + tail, sign_bit) + mean_;
                    }

                    GAPP_UNREACHABLE();
                }

                const T z2 = z * z;
                const T f0 = std::exp(T(-0.5) * (x[i] * x[i] - z2));
                const T f1 = std::exp(T(-0.5) * (x[i + 1] * x[i + 1] - z2));
                const T u1 = rng::generate_canonical<T>(gen);

                if (f1 + u1 * (f0 - f1) < T(1.0))
                {
                    return stddev_ * detail::set_sign_bit(z, sign_bit) + mean_;
                }
            }

            GAPP_UNREACHABLE();
        }

        constexpr void reset() noexcept {}

        constexpr result_type mean() const noexcept { return mean_; }
        constexpr result_type stddev() const noexcept { return stddev_; }

        static constexpr result_type min() noexcept { return -std::numeric_limits<T>::infinity(); }
        static constexpr result_type max() noexcept { return std::numeric_limits<T>::infinity(); }

        friend constexpr bool operator==(const normal_distribution&, const normal_distribution&) = default;

    private:
        T mean_ = 0.0;
        T stddev_ = 1.0;

        static constexpr std::array<T, 129> x =
        {
            3.7130862467425505, 3.44261985589900021, 3.22308498458114157, 3.08322885821686832, 2.97869625264778026,
            2.89434400702152894, 2.82312535054891045, 2.76116937238717686, 2.70611357312181955, 2.65640641126135968,
            2.61097224843184739, 2.56903362592493778, 2.53000967238882746, 2.49345452209537211, 2.45901817741183049,
            2.42642064553374981, 2.39543427801106246, 2.36587137011763859, 2.33757524133923678, 2.31041368369876299,
            2.28427405967747177, 2.25905957386919853, 2.23468639559097948, 2.21108140887870341, 2.18818043207604918,
            2.16592679374892194, 2.14427018236039535, 2.12316570867397658, 2.1025731351892385,  2.08245623799201685,
            2.06278227450830842, 2.04352153665506764, 2.02464697337738553, 2.00613386996347209, 1.98795957412761992,
            1.97010326085432652, 1.95254572955355665, 1.93526922829662285, 1.91825730086450985, 1.90149465310515109,
            1.88496703570775903, 1.86866114099448866, 1.85256451172809111, 1.83666546025844601, 1.82095299659612553,
            1.80541676421922848, 1.79004698259985862, 1.77483439558606948, 1.75977022489959345, 1.74484612811380035,
            1.73005416056373051, 1.71538674071366759, 1.70083661856991686, 1.68639684677916812, 1.67206075409760091,
            1.65782192095402414, 1.64367415686286855, 1.62961147947063467, 1.61562809504316096, 1.60171838022137814,
            1.58787686489057611, 1.57409821602300082, 1.56037722236616894, 1.5467087798599104, 1.53308787767404331,
            1.51950958476594011, 1.50596903686320327, 1.49246142378135405, 1.47898197698992417, 1.46552595734271085,
            1.45208864288922457, 1.43866531668456354, 1.4252512545140601, 1.4118417124470577, 1.39843191413100532,
            1.38501703773265183, 1.37159220242734259, 1.35815245433014353, 1.34469275175354697, 1.33120794966562728,
            1.31769278320941408, 1.30414185012861683, 1.29054959192619645, 1.27691027356015563, 1.26321796145462106,
            1.24946649957306821, 1.23564948326336266, 1.22176023053999638, 1.20779175041594966, 1.19373670783312869,
            1.17958738466398816, 1.16533563616475244, 1.15097284214886741, 1.13648985201316077, 1.12187692258254224,
            1.10712364753403603, 1.09221887690727737, 1.0771506248928957, 1.06190596369482426, 1.04647090076404536,
            1.03083023606819557, 1.01496739525133051, 0.998864233492983589, 0.982500803515429011, 0.965855079401149896,
            0.948902625511306441, 0.931616196615150827, 0.913965251023032277, 0.895915352580937685, 0.877427429112923374,
            0.858456843193813213, 0.838952214297577381, 0.818853906700357292, 0.798092060644056911, 0.776583987894759908,
            0.754230664454055622, 0.730911910642488838, 0.706479611335436464, 0.680747918669154628, 0.653478638739975248,
            0.624358597336050702, 0.592962942471448318, 0.558692178408185192, 0.520656038762060569, 0.477437837296689815,
            0.426547986355423514, 0.36287143109703196, 0.272320864813964669, 0.0
        };
    };


    /** Generates random floating point number from a normal distribution. */
    template<std::floating_point T>
    class normal_distribution_polar
    {
    public:
        static_assert(sizeof(T) <= sizeof(uint64_t));

        using result_type = T;

        constexpr normal_distribution_polar() noexcept = default;

        constexpr normal_distribution_polar(T mean, T stddev) noexcept :
            mean_(mean), stddev_(stddev)
        {
            GAPP_ASSERT(stddev_ > 0.0);
        }

        template<rng::prng64 G>
        constexpr result_type operator()(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
        {
            // Marsaglia polar method

            if (saved_.has_value())
            {
                return *std::exchange(saved_, {}) * stddev_ + mean_;
            }

            for (;;)
            {
                const T x = rng::generate_canonical<T>(gen) - T(0.5);
                const T y = rng::generate_canonical<T>(gen) - T(0.5);
                const T r = x * x + y * y;

                if (r > T(0.25) || r == T(0.0))
                    continue;

                const T scale = std::sqrt(T(-2.0) * std::log(T(4.0) * r) / r);

                const T v1 = x * scale;
                const T v2 = y * scale;

                saved_ = v2;
                return v1 * stddev_ + mean_;
            }

            GAPP_UNREACHABLE();
        }

        constexpr void reset() noexcept { saved_.reset(); }

        constexpr result_type mean() const noexcept { return mean_; }
        constexpr result_type stddev() const noexcept { return stddev_; }

        static constexpr result_type min() noexcept { return -std::numeric_limits<T>::infinity(); }
        static constexpr result_type max() noexcept { return std::numeric_limits<T>::infinity(); }

        friend constexpr bool operator==(const normal_distribution_polar&, const normal_distribution_polar&) = default;

    private:
        T mean_ = 0.0;
        T stddev_ = 1.0;
        std::optional<T> saved_ = {};
    };


    /** Generates random integers from a poisson distribution. Intended for small mean values (<= 16). */
    template<std::integral T>
    class small_poisson_distribution
    {
    public:
        static_assert(sizeof(T) <= sizeof(uint64_t));

        using result_type = T;

        explicit constexpr small_poisson_distribution(double mean) noexcept :
            mean_(mean), mean_exp_(std::exp(-mean))
        {
            GAPP_ASSERT(mean > 0.0);
        }

        template<rng::prng64 G>
        constexpr result_type operator()(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
        {
            T k = 0;
            double pdf = mean_exp_;
            double cdf = rng::generate_canonical<double>(gen);

            while (cdf > pdf)
            {
                cdf = cdf - pdf;
                pdf = pdf * mean_ / ++k;
            }

            return k;
        }

        constexpr void reset() noexcept {}

        constexpr result_type min() const noexcept { return 0; }
        constexpr result_type max() const noexcept { return std::numeric_limits<T>::max(); }

        friend constexpr bool operator==(const small_poisson_distribution&, const small_poisson_distribution&) = default;

    private:
        double mean_;
        double mean_exp_;
    };


    /** Generates random integer from a symmetric binomial distribution (p = 0.5). */
    template<std::integral T>
    class symmetric_binomial_distribution
    {
    public:
        static_assert(sizeof(T) <= sizeof(uint64_t));

        using result_type = T;

        explicit constexpr symmetric_binomial_distribution(T n) noexcept :
            n_(n)
        {
            GAPP_ASSERT(n >= 0);

            if (n_ > 256)
            {
                const double mean = 0.5 * n_;
                const double sdev = 0.5 * std::sqrt(n_);
                norm_ = rng::normal_distribution<double>(mean, sdev);
            }
        }

        template<rng::prng64 G>
        constexpr result_type operator()(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
        {
            if (n_ <= 256)
            {
                const uint64_t div = n_ / 64;
                const uint64_t rem = n_ % 64;

                result_type k = 0;
                for (uint64_t i = 0; i < div; i++)
                {
                    k += result_type(std::popcount(gen()));
                }

                if (!rem) return k;

                k += result_type(std::popcount(gen() & detail::mask_left_n<uint64_t>(rem)));
                return k;
            }

            for (;;)
            {
                const double k = norm_(gen) + 0.5;
                if (k < 0.0 || k > n_)
                    continue;

                return result_type(k);
            }

            GAPP_UNREACHABLE();
        }

        constexpr void reset() noexcept {}

        constexpr result_type min() const noexcept { return 0; }
        constexpr result_type max() const noexcept { return n_; }

        friend constexpr bool operator==(const symmetric_binomial_distribution& lhs, const symmetric_binomial_distribution& rhs) noexcept
        {
            return lhs.n_ == rhs.n_;
        }

    private:
        rng::normal_distribution<double> norm_;
        uint64_t n_;
    };


    /** Generates random integers from a binomial distribution. */
    template<std::integral T>
    class binomial_distribution
    {
    public:
        static_assert(sizeof(T) <= sizeof(uint64_t));

        using result_type = T;

        constexpr binomial_distribution(T n, double p) noexcept :
            n_(n), p_(p), mirror_(p > 0.5)
        {
            GAPP_ASSERT(n >= 0);
            GAPP_ASSERT(0.0 <= p && p <= 1.0);

            if (p == 0.5)
            {
                symm_ = symmetric_binomial_distribution<T>(n);
                return;
            }

            if (mirror_)
            {
                p_ = 1.0 - p_;
            }

            invert_ = n_ * p_ <= 16.0;
            if (invert_)
            {
                binv_.qn = std::pow(1.0 - p_, n_);
                binv_.pdq = p_ / (1.0 - p_);
                binv_.pdqn = binv_.pdq * (n_ + 1.0);
                return;
            }

            btrs_.spq = std::sqrt(n_ * p_ * (1.0 - p_));
            btrs_.lpq = std::log(p_ / (1.0 - p_));
            btrs_.b = 1.15 + 2.53 * btrs_.spq;
            btrs_.a = -0.0873 + 0.0248 * btrs_.b + 0.01 * p_;
            btrs_.c = n_ * p_ + 0.5;
            btrs_.alpha = (2.83 + 5.1 / btrs_.b) * btrs_.spq;
            btrs_.vr = 0.92 - 4.2 / btrs_.b;
            btrs_.m = T((n_ + 1.0) * p_);
            btrs_.h = log_factorial(btrs_.m) + log_factorial(n_ - btrs_.m);
        }

        template<rng::prng64 G>
        constexpr result_type operator()(G& gen) noexcept(std::is_nothrow_invocable_v<G>)
        {
            if (p_ == 0.5) return symm_(gen);

            // Inverse transform algorithm, based on:
            //  Kachitvichyanukul, Voratas, and Bruce W. Schmeiser. "Binomial random variate generation."
            //  Communications of the ACM 31, no. 2 (1988): 216-222.

            if (invert_)
            {
                T k = 0;
                double pdf = binv_.qn;
                double cdf = rng::generate_canonical<double>(gen);

                while (cdf > pdf)
                {
                    cdf = cdf - pdf;
                    pdf = pdf * (binv_.pdqn / ++k - binv_.pdq);
                }

                return mirror_ ? n_ - k : k;
            }

            // BTRS algorithm, based on:
            //  Hörmann, Wolfgang. "The generation of binomial random variates."
            //  Journal of statistical computation and simulation 46, no. 1-2 (1993): 101-110.

            for (;;)
            {
                const double u = rng::generate_canonical<double>(gen) - 0.5;
                const double v = rng::generate_canonical<double>(gen);

                const double us = 0.5 - std::abs(u);
                const double kf = (2.0 * btrs_.a / us + btrs_.b) * u + btrs_.c;

                if (kf < 0.0 || kf >= n_ + 1.0)
                {
                    continue;
                }

                const result_type k = kf;

                if (us >= 0.07 && v <= btrs_.vr)
                {
                    return mirror_ ? n_ - k : k;
                }

                const double v2 = std::log(v * btrs_.alpha / (btrs_.a / (us * us) + btrs_.b));
                const double t = btrs_.h - log_factorial(k) - log_factorial(n_ - k) + (k - btrs_.m) * btrs_.lpq;

                if (v2 <= t)
                {
                    return mirror_ ? n_ - k : k;
                }
            }

            GAPP_UNREACHABLE();
        }

        constexpr void reset() noexcept {}

        constexpr result_type min() const noexcept { return 0; }
        constexpr result_type max() const noexcept { return n_; }

        friend constexpr bool operator==(const binomial_distribution& lhs, const binomial_distribution& rhs) noexcept
        {
            return lhs.n_ == rhs.n_ && 
                   lhs.p_ == rhs.p_ &&
                   lhs.mirror_ == rhs.mirror_;
        }

    private:
        uint64_t n_;
        double p_;

        struct btrs_params
        {
            double spq;
            double lpq;
            double a, b, c;
            double alpha;
            double vr;
            double h;
            T m;
        };

        struct binv_params
        {
            double qn;
            double pdq;
            double pdqn;
        };

        union
        {
            btrs_params btrs_;
            binv_params binv_;
            symmetric_binomial_distribution<T> symm_;
        };

        bool mirror_ = false;
        bool invert_ = false;

        static constexpr std::array<double, 10> log_factorial_table =
        {
            0.0,
            0.0,
            0.6931471805599453,
            1.791759469228055,
            3.1780538303479458,
            4.787491742782046,
            6.579251212010101,
            8.525161361065415,
            10.60460290274525,
            12.801827480081469
        };

        static constexpr double log_factorial(T k) noexcept
        {
            if (k < 10) return log_factorial_table[k];

            constexpr double log_sqrt_2pi = 0.9189385332046727;
            
            const double k1 = k + 1.0;
            const double k1_inv = 1.0 / (k + 1.0);
            const double k1_inv_sq = k1_inv * k1_inv;

            return log_sqrt_2pi + (k + 0.5) * std::log(k1) - k1 + 
                (1.0 / 12.0 - (1.0 / 360.0 - (1.0 / 1260.0 * k1_inv_sq)) * k1_inv_sq) * k1_inv;
        }
    };

} // namespace gapp::rng

#endif // !GAPP_UTILITY_DISTRIBUTION_HPP
