#include "type_def.h"
#include <vector>
#include <algorithm>
#include <iostream>

namespace {

    struct Item {
        int id;
        std::vector<int> tidset;
    };


    // =====================================================
    // intersection (Eclat核心)
    // =====================================================
    std::vector<int> intersect(const std::vector<int>& a,
        const std::vector<int>& b)
    {
        std::vector<int> res;
        std::set_intersection(
            a.begin(), a.end(),
            b.begin(), b.end(),
            std::back_inserter(res)
        );
        return res;
    }

    // =====================================================
    // maximal checking（局部判定版，避免全局扫描）
    // =====================================================
    bool is_subset(const std::vector<int>& a,
        const std::vector<int>& b)
    {
        return std::includes(b.begin(), b.end(), a.begin(), a.end());
    }

    // =====================================================
    // Eclat DFS
    // =====================================================
    void eclat(
        const std::vector<Item>& items,
        int min_sup,
        size_t start,
        std::vector<int> prefix_items,
        std::vector<int> prefix_tidset,
        std::vector<Pattern>& results
    )
    {
        bool has_extension = false;

        for (size_t i = start; i < items.size(); ++i) {

            // 1. compute intersection (tidset)
            std::vector<int> new_tidset =
                intersect(prefix_tidset, items[i].tidset);

            // 2. support pruning
            if (new_tidset.size() < (size_t)min_sup)
                continue;

            has_extension = true;

            // 3. extend pattern
            auto new_items = prefix_items;
            new_items.push_back(items[i].id);

            // 4. recursive DFS
            eclat(items,
                min_sup,
                i + 1,
                new_items,
                new_tidset,
                results);
        }

        // 5. maximal pattern condition
        if (!has_extension) {
            results.push_back({
                prefix_items,
                prefix_tidset,
                (int)prefix_tidset.size()
                });
        }
    }

} // namespace

// =====================================================
// main entry
// =====================================================
void pattern_mining(const KMBCPInstance& instance,
    std::vector<Pattern>& mp_res)
{
    std::vector<Item> items;
    items.reserve(instance.n);
    Pattern test;

    // -------------------------------------------------
    // build vertical tidset
    // -------------------------------------------------
    for (int v = 0; v < instance.n; ++v) {

        std::vector<int> tidset;

        for (int c = 0; c < instance.m; ++c) {
            if (instance.edges[v][instance.n + c]) {
                tidset.push_back(c);
            }
        }

        if (tidset.size() >= (size_t)Paras_gb.min_sup) {
            items.push_back({ v, tidset });
        }
    }

    // sort tidset (required for intersection correctness)
    for (auto& it : items) {
        std::sort(it.tidset.begin(), it.tidset.end());
    }

    mp_res.clear();

    // -------------------------------------------------
    // start Eclat DFS
    // -------------------------------------------------
    for (size_t i = 0; i < items.size(); ++i) {

        eclat(
            items,
            Paras_gb.min_sup,
            i + 1,
            { items[i].id },
            items[i].tidset,
            mp_res
        );
    }

    // -------------------------------------------------
    // output
    // -------------------------------------------------
    std::cout << "Maximal Frequent Patterns:\n";

    for (auto& p : mp_res) {
        std::cout << "{ ";
        for (int x : p.items)
            std::cout << (x + 1) << " ";
        std::cout << "} sup=" << p.support << "\n";
    }
}