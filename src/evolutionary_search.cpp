#include"type_def.h"
#include"common_func.h"
#include<numeric>
#include <unordered_set>
#include <unordered_map>



// 可选：给模式打一个综合分数（支持度、长度都重要）
static inline double pattern_score(const Pattern& p, double alpha = 1.0, double beta = 0.5) {
    // score = support^alpha * |items|^beta
    return pow((double)std::max(1, p.support), alpha) * pow((double)std::max(1, (int)p.items.size()), beta);
}

// 将 cluster 编码为字符串，方便哈希去重
std::string encodeCluster(const std::vector<int>& cluster) {
    std::string key;
    key.reserve(cluster.size() * 2);
    for (int c : cluster) {
        key.push_back(char('A' + (c % 26))); // 简单编码
        key.push_back(',');
    }
    return key;
}


/*
* 将 cluster 转为标准化表示，用于去重
* 这方式真挺好的， GPT厉害啊。
*/
std::string canonicalCluster(const std::vector<int>& cluster, int K) {
    std::unordered_map<int, std::vector<int>> clustersMap;
    for (int i = 0; i < cluster.size(); ++i)
        clustersMap[cluster[i]].push_back(i);

    std::vector<std::vector<int>> clusters;
    for (auto& [c, vec] : clustersMap) {
        sort(vec.begin(), vec.end());
        clusters.push_back(vec);
    }
    sort(clusters.begin(), clusters.end());

    std::stringstream ss;
    for (auto& vec : clusters) {
        for (int v : vec) ss << v << ",";
        ss << ";";
    }
    return ss.str();
}


// 单点移动邻域的快速下降局部搜索
void descent_local_search(Individual& ind, const KMBCPInstance& G, std::mt19937& gen)
{
    Solution sol(G.n, G.m, G.K);
    std::vector<int> u_nodes(G.n);
    bool is_improved = true;
    Help_struc_move struc_move(G.n, G.K);

    iota(u_nodes.begin(), u_nodes.end(), 0);   // u_nodes = {0, 1, ..., n - 1}
    shuffle(u_nodes.begin(), u_nodes.end(), gen);
    // 输出结果
    //for (int x : u_nodes) std::cout << x << " ";
    //std::cout << std::endl;
    //getchar();

    sol.build_solution_data(G, ind.cluster);
   
    while (is_improved) {
        is_improved = false;
        for (int node: u_nodes) {            
            int currK = sol.serverCluster[node];
            for (int k = 0; k < G.K; ++k) {
                if (k == currK || sol.sizeCluster[currK] == 1)
                    continue;
                int delta = computeMoveDelta_improved(G, sol, struc_move, node, currK, k);
                //int delta = computeMoveDelta(G, sol, node, currK, k);

                if (delta < 0) {
                    is_improved = true;
                    applyMove(G, sol, struc_move, node, k, delta);
                    break;
                }
            }
            if (is_improved) {
                break;
            }
        }
    }
    for (int i = 0; i < G.n; i++) {
        ind.cluster[i] = sol.serverCluster[i];
    }
    ind.fitness = sol.cost;
    //getchar();
}

// 基于模式+随机+局部搜索的种群初始化
void initialize_population_with_patterns_and_local_search(
    std::vector<Individual>& population,
    int pop_size, const KMBCPInstance& G,
    const std::vector<Pattern>& patterns,
    std::mt19937& gen, std::unordered_set<std::string>& seen)
{
    population.clear();
    std::uniform_int_distribution<int> rndK(0, G.K - 1);

    // 模式排序（高分优先）
    std::vector<int> order(patterns.size());
    iota(order.begin(), order.end(), 0);
    //stable_sort稳定排序，通过lambda比较器
    stable_sort(order.begin(), order.end(), [&](int a, int b) {
        double sa = pattern_score(patterns[a]);
        double sb = pattern_score(patterns[b]);
        if (sa != sb) return sa > sb;
        if (patterns[a].support != patterns[b].support) return patterns[a].support > patterns[b].support;
        return patterns[a].items.size() > patterns[b].items.size();
        });


    // 1. 利用模式生成初始个体
    int pidx = 0;
    while (population.size() < pop_size && pidx < pop_size * 2) { // 尝试生成 twice pop_size 次
        Individual ind(G.n);
        std::vector<int> cluster_load(G.K, 0);

        // 随机打乱模式顺序
        std::vector<int> shuffled = order;
        shuffle(shuffled.begin(), shuffled.end(), gen);                 //注释该行， 不打乱已排好序的模式， order. 否则，模式随机打乱， 也就是order之前的排序不起作用了。

        for (int id : shuffled) {
            const auto& pat = patterns[id];
            int overlap = 0;
            for (int u : pat.items)
                if (u >= 0 && u < G.n && ind.cluster[u] != -1) overlap++;
            if (overlap > 0) continue;

            int minLoad = *min_element(cluster_load.begin(), cluster_load.end());
            std::vector<int> candidates;
            for (int c = 0; c < G.K; ++c)
                if (cluster_load[c] == minLoad) candidates.push_back(c);
            std::uniform_int_distribution<int> pick(0, static_cast<int>(candidates.size()) - 1);
            int bestC = candidates[pick(gen)];

            for (int u : pat.items) ind.cluster[u] = bestC;
            cluster_load[bestC] += (int)pat.items.size();
        }

        // 填充剩余顶点
        for (int u = 0; u < G.n; ++u)
            if (ind.cluster[u] == -1) ind.cluster[u] = rndK(gen);

        // 局部搜索提升个体
        descent_local_search(ind, G, gen);

        // 去重
        std::string key = canonicalCluster(ind.cluster, G.K);
        if (!seen.count(key)) {
            seen.insert(key);
            if (ind.fitness < Vars_gb.Best_sol.cost) {
                Vars_gb.Best_sol.cost = ind.fitness;
                Vars_gb.Best_sol.serverCluster = ind.cluster;
                auto now = std::chrono::steady_clock::now();
                double elapsed_time = std::chrono::duration<double>(now - Vars_gb.start_time).count();
                Vars_gb.first_time = elapsed_time;
            }
            population.push_back(std::move(ind));
        }

        ++pidx;
    }

    pidx = 0;
    // 2. 补充随机个体
    while (population.size() < pop_size && pidx < pop_size * 2) {
        Individual ind(G.n);
        for (int u = 0; u < G.n; ++u)
            ind.cluster[u] = rndK(gen);

        descent_local_search(ind, G, gen); // 提升个体

        std::string key = canonicalCluster(ind.cluster, G.K);
        if (!seen.count(key)) {
            seen.insert(key);
            if (ind.fitness < Vars_gb.Best_sol.cost) {
                Vars_gb.Best_sol.cost = ind.fitness;
                Vars_gb.Best_sol.serverCluster = ind.cluster;
                auto now = std::chrono::steady_clock::now();
                double elapsed_time = std::chrono::duration<double>(now - Vars_gb.start_time).count();
                Vars_gb.first_time = elapsed_time;
            }
            population.push_back(std::move(ind));
        }
        pidx++;
    }

    std::cout << "Final population size: " << population.size() << std::endl;
}

// 从父代提取“父代中真正被使用”的模式（即该模式内所有 u 的簇号一致）
// 返回：向量< (items, cluster_id) >
static std::vector<std::pair<std::vector<int>, int>> extract_used_patterns(
    const Individual& P, const std::vector<Pattern>& patterns, int U_size)
{
    std::vector<std::pair<std::vector<int>, int>> used;
    for (const auto& pat : patterns) {
        int c = -1; bool ok = true;
        for (int u : pat.items) {
            if (u < 0 || u >= U_size) continue;
            if (P.cluster[u] == -1) { ok = false; break; }
            if (c == -1) c = P.cluster[u];
            else if (P.cluster[u] != c) { ok = false; break; }
        }
        if (ok && c != -1) used.push_back({ pat.items, c });
    }
    return used;
}

Individual crossover_pattern_preserving(const std::vector<Individual>& pop, const KMBCPInstance& ins,
    const std::vector<Pattern>& patterns, std::mt19937& gen)
{
    std::uniform_int_distribution<int> rndK(0, ins.K - 1);
    std::uniform_int_distribution<int> rndP(0, static_cast<int>(pop.size()) - 1);
    Individual P1(ins.n), P2(ins.n);
    Individual C(ins.n);
    C.cluster.assign(ins.n, -1);
    C.fitness = 0;
    if (pop.size() < 2)
        return pop[0];
    int p1 = rndP(gen);
    int p2 = rndP(gen);
    while (p1 == p2)
        p2 = rndP(gen);
    P1 = pop[p1];
    P2 = pop[p2];

    // 1) 收集两个父代中“成块使用”的模式，并按分数排序
    auto used1 = extract_used_patterns(P1, patterns, ins.n);
    auto used2 = extract_used_patterns(P2, patterns, ins.n);

    // 合并（父1优先可减少冲突；也可按模式分数排序）
    std::vector<std::pair<std::vector<int>, int>> blocks;
    blocks.reserve(used1.size() + used2.size());
    blocks.insert(blocks.end(), used1.begin(), used1.end());
    blocks.insert(blocks.end(), used2.begin(), used2.end());

    // 也可按模式分数排序，让大且高支持度的块先落盘
    auto score_of = [&](const std::vector<int>& items) {
        // 找到对应模式打分（简单起见：长度）
        return (int)items.size();
    };
    stable_sort(blocks.begin(), blocks.end(),
        [&](const auto& A, const auto& B) { return score_of(A.first) > score_of(B.first); });

    // 2) 逐块“冻结”到子代（冲突顶点遵循“先到先得”，或用投票）
    std::vector<char> fixed(ins.n, 0);
    for (auto& blk : blocks) {
        int c = blk.second;
        // 统计该块与当前子代的冲突量，若太大可跳过/部分采用
        int conflict = 0;
        for (int u : blk.first) if (u >= 0 && u < ins.n && fixed[u] && C.cluster[u] != c) conflict++;
        if (conflict > 0) continue; // 简化：发生冲突则跳过该块（也可细化到只放不冲突的顶点）

        for (int u : blk.first) {
            if (u < 0 || u >= ins.n) continue;
            if (!fixed[u]) { C.cluster[u] = c; fixed[u] = 1; }
        }
        /*
        // 遍历块内每个顶点
        for (int u : blk.first) {
            if (u < 0 || u >= ins.n) continue;

            // 仅将未定顶点放入子代
            if (!fixed[u]) {
                C.cluster[u] = c;
                fixed[u] = 1;
            }
            // 如果已定，则直接跳过，不丢弃整块
        }*/

    }

    // 3) 对未定顶点：从两亲本继承（多数投票），仍未定则随机
    for (int u = 0; u < ins.n; ++u) {
        if (fixed[u]) continue;
        int a = P1.cluster[u], b = P2.cluster[u];
        if (a != -1 && b != -1) {
            C.cluster[u] = (a == b ? a : (gen() & 1 ? a : b));
        }
        else if (a != -1) C.cluster[u] = a;
        else if (b != -1) C.cluster[u] = b;
        else C.cluster[u] = rndK(gen);
    }

    return C;
}

/*
* 生成反向解（保留模式结构 + 随机扰动部分点）
*/
Individual generate_opposite(const Individual& C,
    const KMBCPInstance& ins,
    const std::vector<Pattern>& patterns,
    std::mt19937& gen,
    double perturb_rate = 0.1) // 10% 顶点扰动
{
    Individual Opp(ins.n);
    Opp.cluster = C.cluster; // 先复制原始解

    std::uniform_real_distribution<double> prob(0.0, 1.0);
    std::uniform_int_distribution<int> rndK(0, ins.K - 1);

    // 1) 识别模式内顶点
    std::unordered_set<int> pattern_vertices;
    for (const auto& pat : patterns) {
        for (int u : pat.items) pattern_vertices.insert(u);
    }

    // 2) 对非模式顶点进行随机扰动
    for (int u = 0; u < ins.n; ++u) {
        if (pattern_vertices.count(u) == 0) { // 非模式顶点
            if (prob(gen) < perturb_rate) {
                int oldc = Opp.cluster[u];
                int newc;
                do { newc = rndK(gen); } while (newc == oldc);
                Opp.cluster[u] = newc;
            }
        }
    }

    // 3) 对模式顶点，可随机扰动少数点，增强差异
    for (int u : pattern_vertices) {
        if (prob(gen) < perturb_rate * 0.5) { // 模式点扰动概率较小
            int oldc = Opp.cluster[u];
            int newc;
            do { newc = rndK(gen); } while (newc == oldc);
            Opp.cluster[u] = newc;
        }
    }

    return Opp;
}


// 变异操作：随机修改若干点的分配
void mutate(Individual& ind, const KMBCPInstance& G, std::mt19937& gen, double mutation_rate = 0.1) {
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    std::uniform_int_distribution<int> rndK(0, G.K - 1);

    for (int i = 0; i < (int)ind.cluster.size(); i++) {
        if (prob(gen) < mutation_rate) {
            // 重新分配该点
            ind.cluster[i] = rndK(gen);
        }
    }
    // 可选：局部搜索修复，使解更合理
    //fast_descend(ind, G);
}

/*
* 增强模式一致性：模式中的顶点尽量被分配到同一个簇。修复或补全解：对局部分散或未分配的顶点进行补齐。
* 提高解的质量：在聚类或分配问题中，模式完整性往往对应更优的目标函数（如减少冲突、平衡负载等）。
* 避免强制移动：通过 ratio 阈值，只对大部分已经聚集的模式顶点进行补全，减少破坏原有良好结构。
*/
void pattern_completion_repair(Individual& C, const std::vector<Pattern>& patterns, int U_size, double ratio = 0.8) {
    for (const auto& pat : patterns) {
        // 统计该模式在各簇的数量
        std::unordered_map<int, int> cnt;
        int valid = 0;
        for (int u : pat.items) {
            if (u < 0 || u >= U_size) continue;
            if (C.cluster[u] != -1) { cnt[C.cluster[u]]++; valid++; }
        }
        if (valid == 0) continue;
        // 找到占比最高的簇
        int bestC = -1, bestCnt = 0;
        for (auto& kv : cnt) { if (kv.second > bestCnt) { bestCnt = kv.second; bestC = kv.first; } }
        if (bestC == -1) continue;
        if ((double)bestCnt >= ratio * (double)pat.items.size()) {
            // 补齐到 bestC
            for (int u : pat.items) {
                if (u < 0 || u >= U_size) continue;
                C.cluster[u] = bestC;
            }
        }
    }
}



// Hamming距离
int hamming_distance(const Individual& A, const Individual& B) {
    int dist = 0;
    for (int i = 0; i < A.cluster.size(); ++i)
        if (A.cluster[i] != B.cluster[i]) dist++;
    return dist;
}

// 种群多样性
double diversity(const Individual& ind, const std::vector<Individual>& population) {
    double sum = 0;
    for (const auto& other : population)
        sum += hamming_distance(ind, other);
    return sum / population.size();
}

// 更新种群：结合适应度、种群多样性，并保证子代唯一
void update_population_with_diversity_unique(std::vector<Individual>& population, const Individual& child, int K, double alpha, std::unordered_set<std::string>& seen) {
    std::string child_key = canonicalCluster(child.cluster, K);
    // 如果子代已经存在，直接舍弃
    if (seen.count(child_key)) return;

    double worst_score = -1e9;
    int replace_idx = -1;

    // 遍历找最差个体
    for (int i = 0; i < population.size(); ++i) {
        double fit_score = -population[i].fitness;
        double div_score = diversity(population[i], population);
        double score = alpha * fit_score + (1.0 - alpha) * div_score;
        if (score > worst_score) {
            worst_score = score;
            replace_idx = i;
        }
    }

    // 子代综合得分
    double child_score = alpha * (-child.fitness) + (1.0 - alpha) * diversity(child, population);

    if (child_score > worst_score && replace_idx != -1) {
        population[replace_idx] = child;
        seen.insert(child_key);
    }
}


//替换掉种群中最差的个体
void update_population_replacing_worst(std::vector<Individual>& pop, const Individual& child, std::unordered_set<std::string>& seen, int K) {
    if (pop.empty())
        return;
    std::string child_key = canonicalCluster(child.cluster, K);
    // 如果子代已经存在，直接舍弃
    if (seen.count(child_key))
        return;

    // 1. 找到最差个体的索引
    int worst_idx = -1;
    double worst_fit = 1e-9;
    for (int i = 0; i < pop.size(); ++i) {
        if (pop[i].fitness > worst_fit) {   // 最小化问题
            worst_fit = pop[i].fitness;
            worst_idx = i;
        }
    }

    // 2. 替换最差个体
    if (worst_idx != -1)
        pop[worst_idx] = child;
}

/*--------------------进化算法主程序-------------------*/
void evolutionary_search(const KMBCPInstance& ins, const std::vector<Pattern>& mp_res, std::mt19937& gen, std::mt19937_64& gen_64)
{
    std::vector<Individual> pop;
    std::unordered_set<std::string> seen;
    Solution child_sol(ins.n, ins.m, ins.K);
    Solution sb(ins.n, ins.m, ins.K);
    Individual child_ind(ins.n);
    Individual child2(ins.n);

    int count_genes = 0;    
    Vars_gb.Best_sol.cost = INT_MAX;

    initialize_population_with_patterns_and_local_search(pop, Paras_gb.pop_size, ins, mp_res, gen, seen);
    for (int i = 0; i < pop.size(); ++i) {
        std::cout << "i=" << i << " ,cost=" << pop[i].fitness << std::endl;
        /*for (int t = 0; t < ins.n; ++t) {
            std::cout << "t:" << t << " ,k=" << pop[i].cluster[t] << std::endl;
        }*/
    }
    //evolution process
    auto now = std::chrono::steady_clock::now();
    double elapsed_time = std::chrono::duration<double>(now - Vars_gb.start_time).count();
    //while (count_genes < generations) {        
    while(elapsed_time < Vars_gb.time_limit){
        child_ind = crossover_pattern_preserving(pop, ins, mp_res, gen);
        pattern_completion_repair(child_ind, mp_res, ins.n, 0.80);
        mutate(child_ind, ins, gen, Paras_gb.mutate_ratio);
      
        child_sol.build_solution_data(ins, child_ind.cluster);
        if ((count_genes + 1) % 200 == 0)
            std::cout << "after crossover operator, child_sol.cost=" << child_sol.cost << std::endl;

        sb = tabuSearch(ins, child_sol, Paras_gb.maxIter_ts, Paras_gb.tabuTenure, gen);

        child2.cluster = sb.serverCluster;
        child2.fitness = sb.cost;
        //update_population_with_diversity_unique(pop, child, instance.K, weight_alpha, seen);
        update_population_replacing_worst(pop, child2, seen, ins.K);

        
        /*
        //opposite-based offspring
        Individual child_opp = generate_opposite(child_ind, ins, mp_res, gen);
        child_sol.build_solution_data(ins, child_opp.cluster);
        if ((count_genes + 1) % 200 == 0)
            std::cout << "after crossover operator, child_opp_sol.cost=" << child_sol.cost << std::endl;

        sb = tabuSearch(ins, child_sol, Paras_gb.maxIter_ts, Paras_gb.tabuTenure, gen);
        child2.cluster = sb.serverCluster;
        child2.fitness = sb.cost;
        //update_population_with_diversity_unique(pop, child, instance.K, weight_alpha, seen);
        update_population_replacing_worst(pop, child2, seen, ins.K);
        */

        count_genes++;     

        now = std::chrono::steady_clock::now();
        elapsed_time = std::chrono::duration<double>(now - Vars_gb.start_time).count();

        std::cout << "count_genes = " << count_genes << ", best_sol.cost = " << Vars_gb.Best_sol.cost << ", first_found_time = " << Vars_gb.first_time <<
            ", elapsed_time = " << elapsed_time << ", after tabu search, cost=" << sb.cost << std::endl;
        std::cout.flush();              //强制刷新缓冲区也是必要的， 特别是在slurm集群上测试时
      
    }
    check_solution(Vars_gb.Best_sol, ins, Vars_gb.Best_sol.cost);
}

