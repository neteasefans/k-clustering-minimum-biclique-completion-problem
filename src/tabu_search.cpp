#include"type_def.h"
#include"common_func.h"
#include <unordered_set>
#include<unordered_map>

struct CandidateMove {
    int i1;
    int i2;
    int k1;
    int k2;
    int type;
};


/*attribute-based tabu search method*/
Solution tabuSearch(const KMBCPInstance& ins, Solution& current_sol, int maxIter, int tabuTenure, std::mt19937& gen)
{
    int iter = 0;
    int non_imp = 0;
    std::vector < std::vector<int >> tabu_list;                   //禁忌表， tabu_list[i][k]代表服务点i，在cluster k中的禁忌期限
    std::vector<CandidateMove> bestMove;
    std::vector<CandidateMove> tabuBestMove;
    std::vector<int> stat_move_types(4, 0);
    std::vector<int> candi_one_moves;
    std::vector<int> candi_swap_moves;
    Help_struc_move struc_move(ins.n, ins.K);
    CandidateMove tmpMove;
    Solution best_sol(ins.n, ins.m, ins.K);
    copy_solution(current_sol, best_sol);
    struc_move.build_strucmove_data(ins, current_sol);
    tabu_list.assign(ins.n, std::vector<int>(ins.K, 0));

    double elapsedTime = getElapsedSeconds(Vars_gb.start_time);

    while (non_imp < maxIter && elapsedTime <= Vars_gb.time_limit) {         //maximum number of non-improvement  consecutive iterations
        int bestDelta = INT_MAX;
        int tabuBestDelta = INT_MAX;
        bestMove.clear();
        tabuBestMove.clear();
        candi_one_moves.clear();
        candi_swap_moves.clear();
        // ====== 1. 枚举所有单点移动 ======         
        for (int i = 0; i < ins.n; ++i) {
            int currK = current_sol.serverCluster[i];
            for (int k = 0; k < ins.K; ++k) {
                if (k == currK || current_sol.sizeCluster[currK] == 1)
                    continue;
                //int delta = computeMoveDelta(ins, current_sol, i, currK, k);
                int delta = computeMoveDelta_improved(ins, current_sol, struc_move, i, currK, k);
                if (tabu_list[i][k] <= iter) {
                    if (delta < bestDelta) {
                        bestDelta = delta;
                        tmpMove = { i, -1, currK, k, 0 };
                        bestMove.clear();
                        bestMove.push_back(tmpMove);
                    }
                    else if (delta == bestDelta) {
                        tmpMove = { i, -1, currK, k, 0};
                        bestMove.push_back(tmpMove);
                    }
                }
                else {
                    if (delta < tabuBestDelta) {
                        tabuBestDelta = delta;
                        tmpMove = { i, -1, currK, k, 0 };
                        tabuBestMove.clear();
                        tabuBestMove.push_back(tmpMove);
                    }
                    else if (delta == tabuBestDelta) {
                        tmpMove = { i, -1, currK, k, 0 };
                        tabuBestMove.push_back(tmpMove);
                    }
                }
            }
        }

        
        // ====== 2. 枚举所有点对交换 ======
        for (int i = 0; i < ins.n; ++i) {
            for (int j = i + 1; j < ins.n; ++j) {
                int ki = current_sol.serverCluster[i];
                int kj = current_sol.serverCluster[j];
                if (ki == kj) continue;
                //int delta = computeSwapDelta(ins, current_sol, i, j);
                int delta = computeSwapDelta_improved(ins, current_sol, struc_move, i, j);
                if (tabu_list[i][kj] <= iter && tabu_list[j][ki] <= iter) {  //not tabu
                    if (delta < bestDelta) {
                        bestDelta = delta;
                        tmpMove = { i, j, ki, kj, 1};
                        bestMove.clear();
                        bestMove.push_back(tmpMove);
                    }
                    else if (delta == bestDelta) {
                        tmpMove = { i, j, ki, kj, 1 };
                        bestMove.push_back(tmpMove);
                    }
                }
                else {
                    if (delta < tabuBestDelta) {
                        tabuBestDelta = delta;
                        tmpMove = { i, j, ki, kj, 1 };
                        tabuBestMove.clear();
                        tabuBestMove.push_back(tmpMove);
                    }
                    else if (delta == tabuBestDelta) {
                        tmpMove = { i, j, ki, kj, 1 };
                        tabuBestMove.push_back(tmpMove);
                    }
                }
            }
        }
        
        // ====== 3. 应用最佳邻域操作 ======
        //aspiration criterion
        if ((tabuBestMove.size() > 0 && tabuBestDelta + current_sol.cost < best_sol.cost && tabuBestDelta < bestDelta) || (tabuBestMove.size() > 0 && bestMove.size() == 0)) {
            std::uniform_int_distribution<int> rndM(0, static_cast<int>(tabuBestMove.size()) - 1);
            int rx = rndM(gen);

            int min_edges = INT_MAX;
            for (int t = 0; t < tabuBestMove.size(); t++) {
                if (tabuBestMove[t].type == 0) {
                    int i = tabuBestMove[t].i1;
                    int fromK = tabuBestMove[t].k1;
                    int toK = tabuBestMove[t].k2;
                    int res = -1 * struc_move.servers_eq_one_intersect_clients[i][fromK].size() - struc_move.servers_ge_one_intersect_clients[i][fromK].size() +
                        struc_move.servers_eq_one_intersect_clients[i][toK].size() + struc_move.servers_ge_one_intersect_clients[i][toK].size();
                    if (res < min_edges) {
                        min_edges = res;
                        candi_one_moves.clear();
                        candi_one_moves.push_back(t);
                    }
                    else if(res == min_edges){
                        candi_one_moves.push_back(t);
                    }
                }
                else {
                    int i = tabuBestMove[t].i1;
                    int i2 = tabuBestMove[t].i2;
                    int k1 = tabuBestMove[t].k1;
                    int k2 = tabuBestMove[t].k2;
                    int res = -1 * struc_move.servers_eq_one_intersect_clients[i][k1].size() - struc_move.servers_ge_one_intersect_clients[i][k1].size() +
                        struc_move.servers_eq_one_intersect_clients[i][k2].size() + struc_move.servers_ge_one_intersect_clients[i][k2].size() -
                        struc_move.servers_eq_one_intersect_clients[i2][k2].size() - struc_move.servers_ge_one_intersect_clients[i2][k2].size() +
                        struc_move.servers_eq_one_intersect_clients[i2][k1].size() + struc_move.servers_ge_one_intersect_clients[i2][k1].size();
                    if (res < min_edges) {
                        min_edges = res;
                        candi_swap_moves.clear();
                        candi_swap_moves.push_back(t);
                    }
                    else if (res == min_edges) {
                        candi_swap_moves.push_back(t);
                    }
                }
            }
            if (candi_one_moves.size() > 0)
                rx = rand() % candi_one_moves.size();
            else
                rx = rand() % candi_swap_moves.size();

            if (tabuBestMove[rx].type == 0) {
                applyMove(ins, current_sol, struc_move, tabuBestMove[rx].i1, tabuBestMove[rx].k2, tabuBestDelta);
                tabu_list[tabuBestMove[rx].i1][tabuBestMove[rx].k1] = iter + tabuTenure;
                stat_move_types[0]++;
            }
            else {
                applySwap(ins, current_sol, struc_move, tabuBestMove[rx].i1, tabuBestMove[rx].i2, tabuBestDelta);
                tabu_list[tabuBestMove[rx].i1][tabuBestMove[rx].k1] = iter + tabuTenure;
                tabu_list[tabuBestMove[rx].i2][tabuBestMove[rx].k2] = iter + tabuTenure;
                stat_move_types[2]++;
            }
        }
        else if(bestMove.size() > 0){           
            std::uniform_int_distribution<int> rndM(0, static_cast<int>(bestMove.size()) - 1);
            int rx = rndM(gen);
            int min_edges = INT_MAX;
            for (int t = 0; t < bestMove.size(); t++) {
                if (bestMove[t].type == 0) {
                    int i = bestMove[t].i1;
                    int fromK = bestMove[t].k1;
                    int toK = bestMove[t].k2;
                    int res = -1 * struc_move.servers_eq_one_intersect_clients[i][fromK].size() - struc_move.servers_ge_one_intersect_clients[i][fromK].size() +
                        struc_move.servers_eq_one_intersect_clients[i][toK].size() + struc_move.servers_ge_one_intersect_clients[i][toK].size();
                    if (res < min_edges) {
                        min_edges = res;
                        candi_one_moves.clear();
                        candi_one_moves.push_back(t);
                    }
                    else if (res == min_edges) {
                        candi_one_moves.push_back(t);
                    }
                }
                else {
                    int i = bestMove[t].i1;
                    int i2 = bestMove[t].i2;
                    int k1 = bestMove[t].k1;
                    int k2 = bestMove[t].k2;
                    int res = -1 * struc_move.servers_eq_one_intersect_clients[i][k1].size() - struc_move.servers_ge_one_intersect_clients[i][k1].size() +
                        struc_move.servers_eq_one_intersect_clients[i][k2].size() + struc_move.servers_ge_one_intersect_clients[i][k2].size() -
                        struc_move.servers_eq_one_intersect_clients[i2][k2].size() - struc_move.servers_ge_one_intersect_clients[i2][k2].size() +
                        struc_move.servers_eq_one_intersect_clients[i2][k1].size() + struc_move.servers_ge_one_intersect_clients[i2][k1].size();
                    if (res < min_edges) {
                        min_edges = res;
                        candi_swap_moves.clear();
                        candi_swap_moves.push_back(t);
                    }
                    else if (res == min_edges) {
                        candi_swap_moves.push_back(t);
                    }
                }
            }
            if (candi_one_moves.size() > 0)
                rx = rand() % candi_one_moves.size();
            else
                rx = rand() % candi_swap_moves.size();

            if (bestMove[rx].type == 0) {
                applyMove(ins, current_sol, struc_move, bestMove[rx].i1, bestMove[rx].k2, bestDelta);
                tabu_list[bestMove[rx].i1][bestMove[rx].k1] = iter + tabuTenure;
                stat_move_types[1]++;
            }
            else {
                applySwap(ins, current_sol, struc_move, bestMove[rx].i1, bestMove[rx].i2, bestDelta);
                tabu_list[bestMove[rx].i1][bestMove[rx].k1] = iter + tabuTenure;
                tabu_list[bestMove[rx].i2][bestMove[rx].k2] = iter + tabuTenure;
                stat_move_types[3]++;
            }
        }
        if (current_sol.cost < best_sol.cost) {
            copy_solution(current_sol, best_sol);
            non_imp = 0;
        }
        if (current_sol.cost < Vars_gb.Best_sol.cost) {
            auto now = std::chrono::steady_clock::now();
            double elapsed_time = std::chrono::duration<double>(now - Vars_gb.start_time).count();
            Vars_gb.first_time = elapsed_time;
            copy_solution(current_sol, Vars_gb.Best_sol);
        }
        iter++;
        non_imp++;
        elapsedTime = getElapsedSeconds(Vars_gb.start_time);

        /*if (iter % 1000 == 0)
            std::cout << "iter=" << iter << ", current_sol.cost=" << current_sol.cost << ", best_sol.cost=" << best_sol.cost
            << ",tabu.size=" << tabuBestMove.size() << " ,nontabu_size= " << bestMove.size() << std::endl;
        check_solution(current_sol, ins, current_sol.cost);
        */
        
        //getchar();       
    }
    //std::cout << "move type stat (type 1: tabu one-move), (type 2: one-move), (type 3: tabu swap-move), (type 4: swap):" << std::endl;
    //for (int i = 0; i < stat_move_types.size(); i++)
     //   std::cout << "move type = " << i + 1 << " , count = " << stat_move_types[i] << std::endl;
    return best_sol;
}



// ----------------------------
// Zobrist 哈希
// ----------------------------
struct ZobristHash {
    int n, K;
    std::vector<std::vector<uint64_t>> Z;
    ZobristHash(const int n_, const int K_, std::mt19937_64& gen) : n(n_), K(K_) {
        Z.assign(n, std::vector<uint64_t>(K));
        for (int i = 0; i < n; i++)
            for (int k = 0; k < K; k++)
                Z[i][k] = gen();
    }

    // 全局哈希
    uint64_t hash(const std::vector<int>& cluster_of) const {
        uint64_t H = 0;
        for (int i = 0; i < n; i++)
            H ^= Z[i][cluster_of[i]];           //^ 为异或XOR操作
        return H;
    }

    // 增量更新：单点移动, O(1) time
    uint64_t update(uint64_t H_old, int idx, int old_k, int new_k) const {
        return H_old ^ Z[idx][old_k] ^ Z[idx][new_k];
    }

    // 多点移动增量更新
    uint64_t update_multiple(uint64_t H_old, const std::vector<std::tuple<int, int, int>>& moves) const {
        uint64_t H = H_old;
        for (auto& [idx, oldk, newk] : moves)
            H ^= Z[idx][oldk] ^ Z[idx][newk];
        return H;
    }
};

// ----------------------------
/* 局部 canonical 化, 消除解的对称性，即(1 1 0 0) 与 （0 0 1 1） 本质是同一个解
* time complexity: O(|C(k1)| + |C(k1)|log(|C(k1)|) + |C(k2)|log(|C(k2)|))，复杂度不低
*/
// ----------------------------
std::vector<int> local_canonical(const std::vector<int>& cluster_of,
    const std::unordered_set<int>& affected_clusters, int K) {
    int n = static_cast<int>(cluster_of.size());
    // 收集受影响簇元素
    std::vector<std::vector<int>> clusters(K);
    for (int i = 0; i < n; i++)
        if (affected_clusters.count(cluster_of[i]))
            clusters[cluster_of[i]].push_back(i);

    // 簇内排序
    for (int k : affected_clusters)
        sort(clusters[k].begin(), clusters[k].end());

    // 生成局部 canonical 映射
    std::vector<int> new_cluster_of = cluster_of;
    std::vector<int> mapping(K, -1);
    int new_id = 0;
    for (int k = 0; k < K; k++) {
        if (affected_clusters.count(k))
            mapping[k] = new_id++;
    }

    for (int i = 0; i < n; i++)
        if (affected_clusters.count(cluster_of[i]))
            new_cluster_of[i] = mapping[cluster_of[i]];

    return new_cluster_of;
}

/*solution-based tabu search method*/
Solution solution_based_tabuSearch(const KMBCPInstance& ins, Solution& current_sol, int maxIter, std::mt19937& gen, std::mt19937_64& gen_64)
{
    int iter = 0;
    int non_imp = 0;
    std::vector<CandidateMove> bestMove;
    std::vector<CandidateMove> tabuBestMove;
    Help_struc_move struc_move(ins.n, ins.K);
    CandidateMove tmpMove;
    Solution best_sol(ins.n, ins.m, ins.K);

    copy_solution(current_sol, best_sol);
    struc_move.build_strucmove_data(ins, current_sol);

    ZobristHash zh(ins.n, ins.K, gen_64);
    uint64_t curr_hash = zh.hash(current_sol.serverCluster);   //the hash value for current solution
    std::unordered_map<uint64_t, int> hash_visited;            /*record the visited solution, by hash value, 第二个元素， 可以记录禁忌期限，本文中未使用禁忌期限*/


    while (non_imp < maxIter) {                                //maximum number of non-improvement  consecutive iterations
        int bestDelta = INT_MAX;
        int tabuBestDelta = INT_MAX;
        bestMove.clear();
        tabuBestMove.clear();

        // ====== 1. 枚举所有单点移动 ======         
        for (int i = 0; i < ins.n; ++i) {
            int currK = current_sol.serverCluster[i];
            for (int k = 0; k < ins.K; ++k) {
                if (k == currK || current_sol.sizeCluster[currK] == 1)
                    continue;
                //int delta = computeMoveDelta(ins, current_sol, i, currK, k);
                int delta = computeMoveDelta_improved(ins, current_sol, struc_move, i, currK, k);
                // 哈希增量更新
                uint64_t new_hash = zh.update(curr_hash, i, currK, k);
                //if (hash_visited.count(new_hash) && hash_visited[new_hash] > iter) {        //tabu
                if (hash_visited.count(new_hash)) {        //tabu
                    if (delta < tabuBestDelta) {
                        tabuBestDelta = delta;
                        tmpMove = { i, -1, currK, k, 0 };
                        tabuBestMove.clear();
                        tabuBestMove.push_back(tmpMove);
                    }
                    else if (delta == tabuBestDelta) {
                        tmpMove = { i, -1, currK, k, 0 };
                        tabuBestMove.push_back(tmpMove);
                    }
                }
                else {
                    if (delta < bestDelta) {
                        bestDelta = delta;
                        tmpMove = { i, -1, currK, k, 0 };
                        bestMove.clear();
                        bestMove.push_back(tmpMove);
                    }
                    else if (delta == bestDelta) {
                        tmpMove = { i, -1, currK, k, 0 };
                        bestMove.push_back(tmpMove);
                    }
                }
            }
        }


        // ====== 3. 应用最佳邻域操作 ======
        //aspiration criterion
        if ((tabuBestMove.size() > 0 && tabuBestDelta + current_sol.cost < best_sol.cost && tabuBestDelta < bestDelta) || (tabuBestMove.size() > 0 && bestMove.size() == 0)) {
            std::uniform_int_distribution<int> rndM(0, static_cast<int>(tabuBestMove.size()) - 1);
            int rx = rndM(gen);

            if (tabuBestMove[rx].type == 0) {
                applyMove(ins, current_sol, struc_move, tabuBestMove[rx].i1, tabuBestMove[rx].k2, tabuBestDelta);
                uint64_t new_hash = zh.update(curr_hash, tabuBestMove[rx].i1, tabuBestMove[rx].k1, tabuBestMove[rx].k2);
                curr_hash = new_hash;
                hash_visited[curr_hash] = iter;
            }
            else {
                applySwap(ins, current_sol, struc_move, tabuBestMove[rx].i1, tabuBestMove[rx].i2, tabuBestDelta);
                std::vector <std::tuple<int, int, int>> moves;
                moves.push_back(std::tuple<int, int, int>(tabuBestMove[rx].i1, tabuBestMove[rx].k1, tabuBestMove[rx].k2));
                moves.push_back(std::tuple<int, int, int>(tabuBestMove[rx].i2, tabuBestMove[rx].k2, tabuBestMove[rx].k1));
                uint64_t new_hash = zh.update_multiple(curr_hash, moves);
                curr_hash = new_hash;
                hash_visited[curr_hash] = iter;
            }
        }
        else if (bestMove.size() > 0) {
            std::uniform_int_distribution<int> rndM(0, static_cast<int>(bestMove.size()) - 1);
            int rx = rndM(gen);

            if (bestMove[rx].type == 0) {
                applyMove(ins, current_sol, struc_move, bestMove[rx].i1, bestMove[rx].k2, bestDelta);
                uint64_t new_hash = zh.update(curr_hash, bestMove[rx].i1, bestMove[rx].k1, bestMove[rx].k2);
                curr_hash = new_hash;
                hash_visited[curr_hash] = iter;
            }
            else {
                applySwap(ins, current_sol, struc_move, bestMove[rx].i1, bestMove[rx].i2, bestDelta);
                std::vector <std::tuple<int, int, int>> moves;
                moves.push_back(std::tuple<int, int, int>(bestMove[rx].i1, bestMove[rx].k1, bestMove[rx].k2));
                moves.push_back(std::tuple<int, int, int>(bestMove[rx].i2, bestMove[rx].k2, bestMove[rx].k1));
                uint64_t new_hash = zh.update_multiple(curr_hash, moves);
                curr_hash = new_hash;
                hash_visited[curr_hash] = iter;
            }
        }
        if (current_sol.cost < best_sol.cost) {
            copy_solution(current_sol, best_sol);
            non_imp = 0;
        }
        if (current_sol.cost < Vars_gb.Best_sol.cost) {
            auto now = std::chrono::steady_clock::now();
            double elapsed_time = std::chrono::duration<double>(now - Vars_gb.start_time).count();
            Vars_gb.first_time = elapsed_time;
            copy_solution(current_sol, Vars_gb.Best_sol);
        }
        iter++;
        non_imp++;
        
        /*
        //if (iter % 1000 == 0)
            std::cout << "iter=" << iter << ", current_sol.cost=" << current_sol.cost << ", best_sol.cost=" << best_sol.cost
            << ",tabu.size=" << tabuBestMove.size() << " ,nontabu_size= " << bestMove.size() << std::endl;
        check_solution(current_sol, ins, current_sol.cost);
        */
    }
    return best_sol;
}
