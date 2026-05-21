#include "test_harness.hpp"

#include <iostream>

int main() {
    std::cout << "\n=== N-Body Unit Tests ===\n\n";

    for ( auto const &t : tests() ) {
        std::cout << "  " << t.name << "... " << std::flush;
        int const fail_before{ g_fail };
        t.fn();
        if ( g_fail == fail_before ) {
            std::cout << "PASS\n";
        }
    }

    std::cout << "\n" << g_pass << " passed, " << g_fail << " failed.\n\n";
    return g_fail > 0 ? 1 : 0;
}
