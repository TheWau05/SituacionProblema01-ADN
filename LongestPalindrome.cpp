#include "Funciones.hpp"
#include <vector>
#include <algorithm>

// See header file for function documentation.
std::string find_longest_palindrome(const std::string& s) {
    if (s.empty()) return "";

    // Transform the string to handle even-length palindromes easily.
    // e.g., "aba" -> "^#a#b#a#$"
    std::string t = "^#";
    for (char c : s) {
        t.push_back(c);
        t.push_back('#');
    }
    t.push_back('$');

    int n = t.size();
    std::vector<int> p(n, 0);
    int center = 0, right = 0;
    int maxLen = 0, centerIndex = 0;

    for (int i = 1; i < n - 1; i++) {
        int mirror = 2 * center - i;
        if (i < right) {
            p[i] = std::min(right - i, p[mirror]);
        }
        // Attempt to expand palindrome centered at i
        while (t[i + 1 + p[i]] == t[i - 1 - p[i]]) {
            p[i]++;
        }
        // If palindrome centered at i expands past `right`,
        // adjust center and right boundary
        if (i + p[i] > right) {
            center = i;
            right = i + p[i];
        }
        // Update the longest palindrome found so far
        if (p[i] > maxLen) {
            maxLen = p[i];
            centerIndex = i;
        }
    }

    int start = (centerIndex - maxLen) / 2;
    return s.substr(start, maxLen);
}
