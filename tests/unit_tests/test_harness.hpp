#pragma once

// Minimal test harness shared across all unit test translation units.
//
// Tests register themselves via the TEST macro using a static initializer.
// inline globals ensure g_pass, g_fail, and the tests() registry are unique
// across the executable (one shared instance across all .cu files).

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <functional>

inline int g_pass{ 0 };
inline int g_fail{ 0 };

struct TestEntry {
    std::string name;
    std::function<void()> fn;
};

// Meyers-singleton pattern: inline function with a static local. All TUs
// see the same vector regardless of which .cu file calls tests().
inline std::vector<TestEntry> &tests() {
    static std::vector<TestEntry> t;
    return t;
}

// Each TEST instantiates a file-scope _reg variable whose initializer
// pushes the test into the global registry at startup.
#define TEST( name ) static void name(); \
    static bool name##_reg = ( tests().push_back( { #name, name } ), true ); \
    static void name()

#define ASSERT_NEAR( actual, expected, tol ) \
    do { \
        double const _a{ (actual) }; \
        double const _e{ (expected) }; \
        double const _t{ (tol) }; \
        if ( std::abs( _a - _e ) > _t ) { \
            std::cerr << "    FAIL: " << #actual << " = " << std::scientific \
                      << std::setprecision( 10 ) << _a << ", expected " << _e \
                      << " +/- " << _t << "\n"; \
            ++g_fail; return; \
        } \
    } while ( 0 )

#define ASSERT_LT( actual, bound ) \
    do { \
        double const _a{ (actual) }; \
        double const _b{ (bound) }; \
        if ( _a >= _b ) { \
            std::cerr << "    FAIL: " << #actual << " = " << std::scientific \
                      << std::setprecision( 10 ) << _a << " >= " << _b << "\n"; \
            ++g_fail; return; \
        } \
    } while ( 0 )

#define ASSERT_TRUE( cond ) \
    do { \
        if ( !(cond) ) { \
            std::cerr << "    FAIL: " << #cond << " is false\n"; \
            ++g_fail; return; \
        } \
    } while ( 0 )
