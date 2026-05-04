#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
    uint64_t x;
    uint64_t y;
    int infinity;
} Point;

typedef struct {
    uint64_t x;
    uint64_t y;
    uint64_t z;
    int infinity;
} JacobianPoint;

typedef struct {
    const char *name;
    uint64_t ecc_prime;
    unsigned __int128 dh_prime;
    int ecc_private_bits;
    int dh_private_bits;
    int iterations;
} Level;

typedef struct {
    double average;
    uint64_t median;
    uint64_t min;
    uint64_t max;
} Result;

static Level levels[] = {
    {"Level 1", 524287ULL, ((unsigned __int128)1 << 61) - 1, 19, 61, 1000},
    {"Level 2", 2147483647ULL, ((unsigned __int128)1 << 89) - 1, 31, 89, 300},
    {"Level 3", 2305843009213693951ULL, ((unsigned __int128)1 << 127) - 1, 61, 127, 50},
};

static uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

static int compare_uint64(const void *a, const void *b) {
    uint64_t left = *(const uint64_t *)a;
    uint64_t right = *(const uint64_t *)b;

    if (left < right) {
        return -1;
    }
    if (left > right) {
        return 1;
    }
    return 0;
}

static uint64_t add_mod(uint64_t a, uint64_t b, uint64_t p) {
    return (uint64_t)(((__uint128_t)a + b) % p);
}

static uint64_t sub_mod(uint64_t a, uint64_t b, uint64_t p) {
    if (a >= b) {
        return a - b;
    }
    return p - (b - a);
}

static uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t p) {
    return (uint64_t)(((__uint128_t)a * b) % p);
}

static uint64_t pow_mod(uint64_t base, uint64_t exponent, uint64_t p) {
    uint64_t result = 1;

    while (exponent > 0) {
        if (exponent & 1) {
            result = mul_mod(result, base, p);
        }

        base = mul_mod(base, base, p);
        exponent >>= 1;
    }

    return result;
}

static unsigned __int128 add_mod128(unsigned __int128 a, unsigned __int128 b,
                                    unsigned __int128 p) {
    return (a + b) % p;
}

static unsigned __int128 mul_mod128(unsigned __int128 a, unsigned __int128 b,
                                    unsigned __int128 p) {
    unsigned __int128 result = 0;

    while (b > 0) {
        if (b & 1) {
            result = add_mod128(result, a, p);
        }

        a = add_mod128(a, a, p);
        b >>= 1;
    }

    return result;
}

static unsigned __int128 pow_mod128(unsigned __int128 base, unsigned __int128 exponent,
                                    unsigned __int128 p) {
    unsigned __int128 result = 1;

    while (exponent > 0) {
        if (exponent & 1) {
            result = mul_mod128(result, base, p);
        }

        base = mul_mod128(base, base, p);
        exponent >>= 1;
    }

    return result;
}

static uint64_t inverse_mod(uint64_t value, uint64_t p) {
    return pow_mod(value, p - 2, p);
}

static int is_on_curve(Point point, uint64_t p) {
    uint64_t left;
    uint64_t right;

    if (point.infinity) {
        return 1;
    }

    left = mul_mod(point.y, point.y, p);
    right = add_mod(mul_mod(mul_mod(point.x, point.x, p), point.x, p), 7, p);

    return left == right;
}

static Point find_generator(uint64_t p) {
    uint64_t x;
    Point point;

    point.x = 0;
    point.y = 0;
    point.infinity = 1;

    for (x = 1; x < 10000; x++) {
        uint64_t right = add_mod(mul_mod(mul_mod(x, x, p), x, p), 7, p);
        uint64_t y = pow_mod(right, (p + 1) / 4, p);

        point.x = x;
        point.y = y;
        point.infinity = 0;

        if (is_on_curve(point, p)) {
            return point;
        }
    }

    return point;
}

/* Jacobian points avoid doing a modular inverse after every point operation. */
static JacobianPoint affine_to_jacobian(Point point) {
    JacobianPoint result;

    result.x = point.x;
    result.y = point.y;
    result.z = point.infinity ? 0 : 1;
    result.infinity = point.infinity;

    return result;
}

static Point jacobian_to_affine(JacobianPoint point, uint64_t p) {
    Point result;
    uint64_t z_inverse;
    uint64_t z2_inverse;
    uint64_t z3_inverse;

    result.x = 0;
    result.y = 0;
    result.infinity = 1;

    if (point.infinity || point.z == 0) {
        return result;
    }

    z_inverse = inverse_mod(point.z, p);
    z2_inverse = mul_mod(z_inverse, z_inverse, p);
    z3_inverse = mul_mod(z2_inverse, z_inverse, p);

    result.x = mul_mod(point.x, z2_inverse, p);
    result.y = mul_mod(point.y, z3_inverse, p);
    result.infinity = 0;

    return result;
}

static JacobianPoint point_double(JacobianPoint point, uint64_t p) {
    JacobianPoint result;
    uint64_t xx;
    uint64_t yy;
    uint64_t yyyy;
    uint64_t s;
    uint64_t m;

    result.x = 0;
    result.y = 0;
    result.z = 0;
    result.infinity = 1;

    if (point.infinity || point.y == 0) {
        return result;
    }

    xx = mul_mod(point.x, point.x, p);
    yy = mul_mod(point.y, point.y, p);
    yyyy = mul_mod(yy, yy, p);
    s = mul_mod(4, mul_mod(point.x, yy, p), p);
    m = mul_mod(3, xx, p);

    result.x = sub_mod(mul_mod(m, m, p), add_mod(s, s, p), p);
    result.y = sub_mod(mul_mod(m, sub_mod(s, result.x, p), p), mul_mod(8, yyyy, p), p);
    result.z = mul_mod(2, mul_mod(point.y, point.z, p), p);
    result.infinity = 0;

    return result;
}

static JacobianPoint point_add(JacobianPoint p1, Point p2, uint64_t p) {
    JacobianPoint result;
    uint64_t z1z1;
    uint64_t u2;
    uint64_t s2;
    uint64_t h;
    uint64_t r;
    uint64_t hh;
    uint64_t hhh;
    uint64_t v;

    if (p1.infinity) {
        return affine_to_jacobian(p2);
    }

    result.x = 0;
    result.y = 0;
    result.z = 0;
    result.infinity = 1;

    if (p2.infinity) {
        return p1;
    }

    z1z1 = mul_mod(p1.z, p1.z, p);
    u2 = mul_mod(p2.x, z1z1, p);
    s2 = mul_mod(p2.y, mul_mod(p1.z, z1z1, p), p);
    h = sub_mod(u2, p1.x, p);
    r = sub_mod(s2, p1.y, p);

    if (h == 0) {
        if (r == 0) {
            return point_double(p1, p);
        }
        return result;
    }

    hh = mul_mod(h, h, p);
    hhh = mul_mod(h, hh, p);
    v = mul_mod(p1.x, hh, p);

    result.x = sub_mod(sub_mod(mul_mod(r, r, p), hhh, p), add_mod(v, v, p), p);
    result.y = sub_mod(mul_mod(r, sub_mod(v, result.x, p), p), mul_mod(p1.y, hhh, p), p);
    result.z = mul_mod(p1.z, h, p);
    result.infinity = 0;

    return result;
}

static Point scalar_multiply(uint64_t scalar, Point point, uint64_t p) {
    JacobianPoint result;
    int bit;

    result.x = 0;
    result.y = 0;
    result.z = 0;
    result.infinity = 1;

    for (bit = 63; bit >= 0; bit--) {
        if (!result.infinity) {
            result = point_double(result, p);
        }

        if ((scalar >> bit) & 1) {
            result = point_add(result, point, p);
        }
    }

    return jacobian_to_affine(result, p);
}

static uint64_t random_bits(int bits) {
    uint64_t value = 0;
    int i;

    for (i = 0; i < 5; i++) {
        value = (value << 15) ^ (uint64_t)(rand() & 0x7fff);
    }

    if (bits < 64) {
        value &= ((1ULL << bits) - 1);
    }

    value |= (1ULL << (bits - 1));
    return value;
}

static unsigned __int128 random_bits128(int bits) {
    unsigned __int128 value = 0;
    int i;

    for (i = 0; i < 9; i++) {
        value = (value << 15) ^ (unsigned __int128)(rand() & 0x7fff);
    }

    if (bits < 128) {
        value &= (((unsigned __int128)1 << bits) - 1);
    }

    value |= ((unsigned __int128)1 << (bits - 1));
    return value;
}

static int ecc_key_exchange(uint64_t p, Point generator, int private_bits) {
    uint64_t alice_private = random_bits(private_bits) % p;
    uint64_t bob_private = random_bits(private_bits) % p;
    Point alice_public = scalar_multiply(alice_private, generator, p);
    Point bob_public = scalar_multiply(bob_private, generator, p);
    Point alice_shared = scalar_multiply(alice_private, bob_public, p);
    Point bob_shared = scalar_multiply(bob_private, alice_public, p);

    if (alice_shared.infinity || bob_shared.infinity) {
        return 0;
    }

    return alice_shared.x == bob_shared.x && alice_shared.y == bob_shared.y;
}

static int regular_dh_key_exchange(unsigned __int128 p, int private_bits) {
    unsigned __int128 generator = 5;
    unsigned __int128 alice_private = random_bits128(private_bits);
    unsigned __int128 bob_private = random_bits128(private_bits);
    unsigned __int128 alice_public = pow_mod128(generator, alice_private, p);
    unsigned __int128 bob_public = pow_mod128(generator, bob_private, p);
    unsigned __int128 alice_shared = pow_mod128(bob_public, alice_private, p);
    unsigned __int128 bob_shared = pow_mod128(alice_public, bob_private, p);

    return alice_shared == bob_shared;
}

static Result make_result(uint64_t times[], int iterations) {
    Result result;
    uint64_t total = 0;
    int i;

    for (i = 0; i < iterations; i++) {
        total += times[i];
    }

    qsort(times, (size_t)iterations, sizeof(uint64_t), compare_uint64);

    result.average = (double)total / (double)iterations;
    result.median = times[iterations / 2];
    result.min = times[0];
    result.max = times[iterations - 1];

    return result;
}

static Result benchmark(int use_ecc, int iterations, uint64_t ecc_p, Point generator,
                        int ecc_private_bits, unsigned __int128 dh_p, int dh_private_bits) {
    Result result;
    uint64_t *times = malloc(sizeof(uint64_t) * (size_t)iterations);
    int i;

    result.average = 0;
    result.median = 0;
    result.min = 0;
    result.max = 0;

    if (times == NULL) {
        return result;
    }

    for (i = 0; i < iterations; i++) {
        uint64_t start;
        uint64_t end;
        int ok;

        if (use_ecc) {
            start = now_ns();
            ok = ecc_key_exchange(ecc_p, generator, ecc_private_bits);
            end = now_ns();
        } else {
            start = now_ns();
            ok = regular_dh_key_exchange(dh_p, dh_private_bits);
            end = now_ns();
        }

        if (!ok) {
            printf("key exchange failed\n");
            free(times);
            return result;
        }

        times[i] = end - start;
    }

    result = make_result(times, iterations);
    free(times);
    return result;
}

static void print_result(const char *name, Result result) {
    printf("%s\n", name);
    printf("  average: %.3f us\n", result.average / 1000.0);
    printf("  median:  %.3f us\n", (double)result.median / 1000.0);
    printf("  min:     %.3f us\n", (double)result.min / 1000.0);
    printf("  max:     %.3f us\n", (double)result.max / 1000.0);
}

static void run_level(Level level, int override_iterations) {
    int iterations = level.iterations;
    Point generator = find_generator(level.ecc_prime);
    Result ecc_result;
    Result dh_result;
    double ratio;

    if (override_iterations > 0) {
        iterations = override_iterations;
    }

    if (generator.infinity) {
        printf("Could not find an ECC generator for %s\n", level.name);
        return;
    }

    printf("%s\n", level.name);
    printf("ECC curve: y^2 = x^3 + 7 mod %llu\n", (unsigned long long)level.ecc_prime);
    printf("ECC generator: (%llu, %llu)\n",
           (unsigned long long)generator.x,
           (unsigned long long)generator.y);
    printf("ECC private key size: %d bits\n", level.ecc_private_bits);
    printf("Regular DH prime size: %d bits\n", level.dh_private_bits);
    printf("Regular DH private key size: %d bits\n", level.dh_private_bits);
    printf("Iterations: %d\n\n", iterations);

    ecc_result = benchmark(1, iterations, level.ecc_prime, generator,
                           level.ecc_private_bits, level.dh_prime, level.dh_private_bits);
    dh_result = benchmark(0, iterations, level.ecc_prime, generator,
                          level.ecc_private_bits, level.dh_prime, level.dh_private_bits);

    print_result("ECC Diffie-Hellman", ecc_result);
    printf("\n");
    print_result("Regular Diffie-Hellman", dh_result);
    printf("\n");

    ratio = ecc_result.average / dh_result.average;
    if (ratio > 1.0) {
        printf("ECC was about %.2f times slower on average.\n", ratio);
    } else {
        printf("ECC was about %.2f times faster on average.\n", 1.0 / ratio);
    }

    printf("This is a self-contained demo benchmark, not production cryptography.\n");
}

int main(int argc, char **argv) {
    int override_iterations = 0;
    size_t i;
    size_t level_count = sizeof(levels) / sizeof(levels[0]);

    srand((unsigned int)time(NULL));

    if (argc > 1) {
        override_iterations = atoi(argv[1]);
        if (override_iterations < 1) {
            printf("Iteration override should be at least 1.\n");
            return 1;
        }
    }

    for (i = 0; i < level_count; i++) {
        run_level(levels[i], override_iterations);

        if (i + 1 < level_count) {
            printf("\n--------------------------------------------------\n\n");
        }
    }

    return 0;
}
