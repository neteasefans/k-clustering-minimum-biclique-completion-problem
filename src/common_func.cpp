#include"type_def.h"
#include<fstream>

/*单点移动Δ计算
* time complexity: O(C(fromK) + C（toK） + d(si)), C(k)为cluster k对应的的客户集合, d(si)为si相邻的客户点集合
* 该增量评估函数， 复杂度较高
*/
int computeMoveDelta(const KMBCPInstance& G, const Solution& sol, int si, int fromK, int toK)
{
    int delta = 0;
    int delta_j = 0;

    //the following is for cluster fromK
    for (int j: sol.clients_in_cluster[fromK]) {
        int deg = sol.clientDegreeInCluster[j][fromK];
        delta_j = 0;
        if (deg > 1) {          
            if(!G.edges[si][j])
                delta_j = -1;
        } 
        else {    //deg == 1
            if (G.edges[si][j])
                delta_j = 1 - sol.sizeCluster[fromK];
            else
                delta_j = -1;
        }            
        delta += delta_j;   
    }

    //the following is for cluster toK
    for (int j: sol.clients_in_cluster[toK]) {
        if (!G.edges[si][j])
            delta += 1;        
    }
    //others， also for toK
    for (int j: G.serverToClients[si]) {
        if (sol.clientDegreeInCluster[j][toK] == 0)
            delta += sol.sizeCluster[toK];      
    }
    //std::cout << "in computeMoveDelta method, detla=" << delta << " , sol.cost=" << sol.cost << ",si=" << si << " ,fromk=" << fromK << " ,toK=" << toK << std::endl;
    return delta;
}

/*单点移动Δ计算, 改进版
* time complexity: O(1)
* 计算效率提升非常明显
* for 15-15-5-0.3.dat, fast incremental evaluation, required time = 0.681 seconds, and not using this strategy, required time = 4.1 sec, 
* 点数越多， 图越稠密， fast incremental evaluation的作用将越发明显。
*/
int computeMoveDelta_improved(const KMBCPInstance& G, const Solution& sol, const Help_struc_move& struc_move, int si, int fromK, int toK)
{
    int delta = -1 * (static_cast<int>(struc_move.servers_ge_one_clients[fromK].size()) - static_cast<int>(struc_move.servers_ge_one_intersect_clients[si][fromK].size())) +
        -1 * (static_cast<int>(struc_move.servers_eq_one_clients[fromK].size()) - static_cast<int>(struc_move.servers_eq_one_intersect_clients[si][fromK].size())) +
        (1 - sol.sizeCluster[fromK]) * (static_cast<int>(struc_move.servers_eq_one_intersect_clients[si][fromK].size()));
    delta += 1 * (static_cast<int>(sol.clients_in_cluster[toK].size()) - static_cast<int>(struc_move.server_adjs_in_clients[si][toK].size()));
    delta += static_cast<int>(struc_move.server_adjs_not_in_clients[si][toK].size()) * sol.sizeCluster[toK];
    
    //std::cout << "in computeMoveDelta method, detla=" << delta << " , sol.cost=" << sol.cost << ",si=" << si << " ,fromk=" << fromK << " ,toK=" << toK << std::endl;
    return delta;
}

/*点对交换Δ计算
* time complexity: O(C(k1)) + C(k2) + d(si) + d(sj)), C(k)为cluster k中客户点集合, d(i)为服务点i的邻接客户集合
* 与单点移动复杂度接近
*/
int computeSwapDelta(const KMBCPInstance& G, const Solution& sol, int si, int sj) 
{
    int delta = 0;
    int delta_j1 = 0;
    int delta_j2 = 0;
    int k1 = sol.serverCluster[si];
    int k2 = sol.serverCluster[sj];

    //the following is for cluster k1
    for (int j: sol.clients_in_cluster[k1]) {
        int deg = sol.clientDegreeInCluster[j][k1];
        delta_j1 = 0;
        delta_j2 = 0;
        if (deg > 1) {            
            if (!G.edges[si][j])
                delta_j1 = -1;
        }
        else {    //deg == 1
            if (G.edges[si][j])
                delta_j1 = 1 - sol.sizeCluster[k1];
            else
                delta_j1 = -1;           
        }

        if (deg > 1) {
            if (!G.edges[sj][j])
                delta_j2 = 1;
        }
        else {
            if (G.edges[si][j]) {
                if (G.edges[sj][j])
                    delta_j2 = sol.sizeCluster[k1] - 1;
            }
            else {
                if (!G.edges[sj][j])
                    delta_j2 = 1;
            }
        }
        delta += delta_j1 + delta_j2;
    }    

    //the following is for cluster k2
    for (int j : sol.clients_in_cluster[k2]) {
        int deg = sol.clientDegreeInCluster[j][k2];
        delta_j1 = 0;
        delta_j2 = 0;
        if (deg > 1) {           
            if(!G.edges[sj][j])
                delta_j1 = -1;
        }
        else {    //deg == 1
            if (G.edges[sj][j])
                delta_j1 = 1 - static_cast<int>(sol.sizeCluster[k2]);
            else
                delta_j1 = -1;
        }

        if (deg > 1) {
            if (!G.edges[si][j])
                delta_j2 = 1;
        }
        else {
            if (G.edges[sj][j]) {
                if (G.edges[si][j])
                    delta_j2 = sol.sizeCluster[k2] - 1;
            }
            else {
                if (!G.edges[si][j])
                    delta_j2 = 1;
            }
        }       
        delta += delta_j1 + delta_j2;
    }

    //others
    for (int j: G.serverToClients[si]) {
        if (sol.clientDegreeInCluster[j][k2] == 0)
            delta += sol.sizeCluster[k2] - 1;
    }
    for (int j: G.serverToClients[sj]) {
        if (sol.clientDegreeInCluster[j][k1] == 0) 
            delta += sol.sizeCluster[k1] - 1;        
    }
    //std::cout << "in computeSwapDelta method, detla=" << delta << " , sol.cost=" << sol.cost << ",si=" << si << " , sj=" << sj << std::endl;
    return delta;
}

/*点对交换Δ计算，改进版
* time complexity: O(eq_one_clients(k1) + eq_one_clients(k2)) 
* 对于稠密图来说， 复杂度应该会比较低， 因为eq_one_clients(k1) + eq_one_clients(k2)会比较少。
*/
int computeSwapDelta_improved(const KMBCPInstance& G, const Solution& sol, const Help_struc_move& struc_move, int si, int sj)
{
    int delta = 0;
    int delta_j = 0;
    int k1 = sol.serverCluster[si];
    int k2 = sol.serverCluster[sj];

    //the following is for cluster k1
    for (int j : struc_move.servers_eq_one_clients[k1]) {
        delta_j = 0;               
        if (G.edges[si][j] && G.edges[sj][j])
            delta_j = sol.sizeCluster[k1] - 1;
        else if (!G.edges[si][j] && (!G.edges[sj][j]))
            delta_j = 1;                    
        delta += delta_j;
    }

    //the following is for cluster k2
    for (int j : struc_move.servers_eq_one_clients[k2]) {
        delta_j = 0;
        if (G.edges[sj][j] && G.edges[si][j])
            delta_j = sol.sizeCluster[k2] - 1;
        else if (!G.edges[sj][j] && (!G.edges[si][j]))
            delta_j = 1;        
        delta += delta_j;
    }

    delta += -1 * (static_cast<int>(struc_move.servers_ge_one_clients[k1].size()) - static_cast<int>(struc_move.servers_ge_one_intersect_clients[si][k1].size())) +
        -1 * (static_cast<int>(struc_move.servers_eq_one_clients[k1].size()) - static_cast<int>(struc_move.servers_eq_one_intersect_clients[si][k1].size())) +
        (1 - sol.sizeCluster[k1]) * (static_cast<int>(struc_move.servers_eq_one_intersect_clients[si][k1].size())) +
        1 * (static_cast<int>(struc_move.servers_ge_one_clients[k1].size()) - static_cast<int>(struc_move.servers_ge_one_intersect_clients[sj][k1].size()));
        
    delta += -1 * (static_cast<int>(struc_move.servers_ge_one_clients[k2].size()) - static_cast<int>(struc_move.servers_ge_one_intersect_clients[sj][k2].size())) +
        -1 * (static_cast<int>(struc_move.servers_eq_one_clients[k2].size()) - static_cast<int>(struc_move.servers_eq_one_intersect_clients[sj][k2].size())) +
        (1 - sol.sizeCluster[k2]) * (static_cast<int>(struc_move.servers_eq_one_intersect_clients[sj][k2].size())) +
        1 * (static_cast<int>(struc_move.servers_ge_one_clients[k2].size()) - static_cast<int>(struc_move.servers_ge_one_intersect_clients[si][k2].size()));
      
    delta += static_cast<int>(struc_move.server_adjs_not_in_clients[si][k2].size()) * (sol.sizeCluster[k2] - 1) +
        static_cast<int>(struc_move.server_adjs_not_in_clients[sj][k1].size()) * (sol.sizeCluster[k1] - 1);
    return delta;
}

void output_res(const std::string outFileName, const std::string input_file, const Solution& best_sol,
    double first_found_time,
    double miningTime,
    int seed)
{
    std::ofstream outfile(outFileName, std::ios::app);
    if (!outfile.is_open()) {
        std::cerr << "Cannot open output file: " << outFileName << std::endl;
        exit(-1);
    }

    // 写入内容
    outfile << input_file << ", best_v = " << best_sol.cost << ", best_time = " << first_found_time <<
        " ,miningTime=" << miningTime << ", seed = " << seed << std::endl;
    for (int k : best_sol.serverCluster)
        outfile << k << " ";
    outfile << std::endl;
    outfile.close();  // 关闭文件
}


/*update struc_move structure after performing a one-move 
* time complexity: O(d(si) ^ 2), d(si)为与服务点si相邻的客户点集合
*/
void update_struc_move_after_onemove(const KMBCPInstance& ins, const Solution& sol, Help_struc_move& struc_move, int si, int fromK, int toK)
{
    for (int j : ins.serverToClients[si]) {
        if (sol.clientDegreeInCluster[j][fromK] == 0) {                 //j is removed from cluster fromK  
            struc_move.servers_eq_one_clients[fromK].erase(std::remove(struc_move.servers_eq_one_clients[fromK].begin(),
                struc_move.servers_eq_one_clients[fromK].end(), j), struc_move.servers_eq_one_clients[fromK].end());

            for (int i: ins.clientToServers[j]) {                
                struc_move.server_adjs_in_clients[i][fromK].erase(std::remove(struc_move.server_adjs_in_clients[i][fromK].begin(),
                    struc_move.server_adjs_in_clients[i][fromK].end(), j), struc_move.server_adjs_in_clients[i][fromK].end());
                struc_move.server_adjs_not_in_clients[i][fromK].push_back(j);

                struc_move.servers_eq_one_intersect_clients[i][fromK].erase(std::remove(struc_move.servers_eq_one_intersect_clients[i][fromK].begin(),
                    struc_move.servers_eq_one_intersect_clients[i][fromK].end(), j), struc_move.servers_eq_one_intersect_clients[i][fromK].end());
            }                     
        }
        if (sol.clientDegreeInCluster[j][fromK] == 1) {
            struc_move.servers_eq_one_clients[fromK].push_back(j);
            struc_move.servers_ge_one_clients[fromK].erase(std::remove(struc_move.servers_ge_one_clients[fromK].begin(),
                struc_move.servers_ge_one_clients[fromK].end(), j), struc_move.servers_ge_one_clients[fromK].end());
                  
            for (int i : ins.clientToServers[j]) {                
                struc_move.servers_eq_one_intersect_clients[i][fromK].push_back(j);
                struc_move.servers_ge_one_intersect_clients[i][fromK].erase(std::remove(struc_move.servers_ge_one_intersect_clients[i][fromK].begin(),
                    struc_move.servers_ge_one_intersect_clients[i][fromK].end(), j), struc_move.servers_ge_one_intersect_clients[i][fromK].end());
            }
        }

        if (sol.clientDegreeInCluster[j][toK] == 1) {

            struc_move.servers_eq_one_clients[toK].push_back(j);
            struc_move.servers_eq_one_intersect_clients[si][toK].push_back(j);

            for (int i : ins.clientToServers[j]) {
                struc_move.server_adjs_not_in_clients[i][toK].erase(std::remove(struc_move.server_adjs_not_in_clients[i][toK].begin(),
                    struc_move.server_adjs_not_in_clients[i][toK].end(), j), struc_move.server_adjs_not_in_clients[i][toK].end());
                struc_move.server_adjs_in_clients[i][toK].push_back(j);
            }
        }
        else if (sol.clientDegreeInCluster[j][toK] == 2) {
            struc_move.servers_eq_one_clients[toK].erase(std::remove(struc_move.servers_eq_one_clients[toK].begin(),
                struc_move.servers_eq_one_clients[toK].end(), j), struc_move.servers_eq_one_clients[toK].end());
            struc_move.servers_ge_one_clients[toK].push_back(j);

            for (int i : ins.clientToServers[j]) {
                struc_move.servers_eq_one_intersect_clients[i][toK].erase(std::remove(struc_move.servers_eq_one_intersect_clients[i][toK].begin(),
                    struc_move.servers_eq_one_intersect_clients[i][toK].end(), j), struc_move.servers_eq_one_intersect_clients[i][toK].end());
                struc_move.servers_ge_one_intersect_clients[i][toK].push_back(j);
            }
        }
    }
}

/*应用单点移动操作
* time complexity: O(d(si) + d(si)^ 2), d(si)为服务点si相邻的客户点集合, d(si)^ 2为update_struc_move 函数的复杂度
*/
void applyMove(const KMBCPInstance& G, Solution& sol, Help_struc_move& struc_move, int si, int toK, int best_delta)
{
    int fromK = sol.serverCluster[si];
    std::vector<int> *vec_fromK = &sol.servers_in_cluster[fromK];
    std::vector<int>* vec_fromK_j = &sol.clients_in_cluster[fromK];

    sol.serverCluster[si] = toK;
    sol.sizeCluster[fromK]--;
    sol.sizeCluster[toK]++;
    for (int j : G.serverToClients[si]) {
        sol.clientDegreeInCluster[j][fromK]--;
        sol.clientDegreeInCluster[j][toK]++;
        if (sol.clientDegreeInCluster[j][fromK] == 0) 
            (*vec_fromK_j).erase(std::remove((*vec_fromK_j).begin(), (*vec_fromK_j).end(), j), (*vec_fromK_j).end());           
              
        if (sol.clientDegreeInCluster[j][toK] == 1) 
            sol.clients_in_cluster[toK].push_back(j);                  
        
    }
    (*vec_fromK).erase(std::remove((*vec_fromK).begin(), (*vec_fromK).end(), si), (*vec_fromK).end());
    sol.servers_in_cluster[toK].push_back(si);
    update_struc_move_after_onemove(G, sol, struc_move, si, fromK, toK);
    sol.cost += best_delta;
}

/*应用点对交换操作*/
void applySwap(const KMBCPInstance& G, Solution& sol, Help_struc_move& struc_move, int si, int sj, int best_delta)
{
    int k1 = sol.serverCluster[si];
    int k2 = sol.serverCluster[sj];

    std::vector<int>* vec_k1 = &sol.servers_in_cluster[k1];
    std::vector<int>* vec_k2 = &sol.servers_in_cluster[k2];
    std::vector<int>* vec_k1_j = &sol.clients_in_cluster[k1];
    std::vector<int>* vec_k2_j = &sol.clients_in_cluster[k2];

    sol.serverCluster[si] = k2;
    sol.serverCluster[sj] = k1;

    for (int j : G.serverToClients[si]) {
        sol.clientDegreeInCluster[j][k1]--;
        sol.clientDegreeInCluster[j][k2]++;
        if (sol.clientDegreeInCluster[j][k1] == 0)
            (*vec_k1_j).erase(std::remove((*vec_k1_j).begin(), (*vec_k1_j).end(), j), (*vec_k1_j).end());

        if (sol.clientDegreeInCluster[j][k2] == 1)
            sol.clients_in_cluster[k2].push_back(j);
    }
    update_struc_move_after_onemove(G, sol, struc_move, si, k1, k2);

    for (int j : G.serverToClients[sj]) {
        sol.clientDegreeInCluster[j][k2]--;
        sol.clientDegreeInCluster[j][k1]++;
        if (sol.clientDegreeInCluster[j][k2] == 0)
            (*vec_k2_j).erase(std::remove((*vec_k2_j).begin(), (*vec_k2_j).end(), j), (*vec_k2_j).end());

        if (sol.clientDegreeInCluster[j][k1] == 1)
            sol.clients_in_cluster[k1].push_back(j);

    }
    update_struc_move_after_onemove(G, sol, struc_move, sj, k2, k1);

    (*vec_k1).erase(std::remove((*vec_k1).begin(), (*vec_k1).end(), si), (*vec_k1).end());
    sol.servers_in_cluster[k2].push_back(si);

    (*vec_k2).erase(std::remove((*vec_k2).begin(), (*vec_k2).end(), sj), (*vec_k2).end());
    sol.servers_in_cluster[k1].push_back(sj);
    sol.cost += best_delta;   
}

/*copy a solution*/
void copy_solution(const Solution& sol_src, Solution& sol_des)
{
    sol_des.serverCluster = sol_src.serverCluster;
    sol_des.sizeCluster = sol_src.sizeCluster;
    sol_des.servers_in_cluster = sol_src.servers_in_cluster;
    sol_des.clientDegreeInCluster = sol_src.clientDegreeInCluster;
    sol_des.clients_in_cluster = sol_src.clients_in_cluster;
    sol_des.cost = sol_src.cost;
}


/*check a move*/
void check_solution(const Solution& sol, const KMBCPInstance& G, int cost)
{
    Solution sol_tmp(G.n, G.m, G.K);
    sol_tmp.serverCluster = sol.serverCluster;
    sol_tmp.sizeCluster.resize(G.K, 0);
    sol_tmp.clientDegreeInCluster.resize(static_cast<size_t>(G.n) + static_cast<size_t>(G.m));
    for (int j = 0; j < G.n + G.m; j++)
        sol_tmp.clientDegreeInCluster[j].resize(G.K, 0);

    for (int i = 0; i < G.n; i++) {
        if (sol.serverCluster[i] < 0 || sol.serverCluster[i] > G.K) {
            std::cout << "in check_solution method, an error is met, sol.serverCluster[k] = " << sol.serverCluster[i] << " ,i= " << i << std::endl;
            getchar();
        }
    }

    for (int i = 0; i < G.n; i++)
        sol_tmp.sizeCluster[sol_tmp.serverCluster[i]]++;

    for (int i = 0; i < G.n; i++) {
        int k = sol_tmp.serverCluster[i];
        for (int j : G.serverToClients[i]) {
            sol_tmp.clientDegreeInCluster[j][k]++;
        }
    }

    sol_tmp.set_cost(G);
    if (sol_tmp.cost != cost) {
        std::cout << "an error is encounterd in check_solution method, cost=" << cost << " , while cost_real=" << sol_tmp.cost << std::endl;
        getchar();
    }
}


/*初始解生成， 随机方式*/
void generateInitialSolutionRandom(const KMBCPInstance& G, Solution& sol, Help_struc_move& struc_move, int seed)
{    
    std::vector<int> tmp_v(G.n);    
    // 用系统时间作为随机种子
    for (int i = 0; i < G.n; i++)
        tmp_v[i] = i;
    // 打乱顺序
    std::shuffle(tmp_v.begin(), tmp_v.end(), std::default_random_engine(seed));   
    for (int i = 0; i < G.n; i++)
        sol.serverCluster[i] = -1;           //初值为-1, 表示该服务点未被分配
    sol.cost = 0;

    // 初始化：先保证每个簇至少有一个节点
    for (int k = 0; k < G.K; k++) {
        int i = tmp_v[k];
        sol.serverCluster[i] = k;
        sol.sizeCluster[k]++;
    }
    for (int i = 0; i < G.n; ++i) {
        if (sol.serverCluster[i] == G.K) {
            int k = rand() % G.K;
            sol.serverCluster[i] = k;           
            sol.sizeCluster[k]++;
        }
    }
    //for testing
    /*sol.serverCluster[0] = 1;
    sol.serverCluster[1] = 1;
    sol.serverCluster[2] = 0;
    sol.serverCluster[3] = 0;
    sol.sizeCluster[0] = 2;
    sol.sizeCluster[1] = 2;*/
    tmp_v = sol.serverCluster;
    sol.build_solution_data(G, tmp_v);
    struc_move.build_strucmove_data(G, sol);
    /*std::cout << " in generateInitialSolutionRandom method ," << std::endl;
    for (int i = 0; i < G.n; i++)
        std::cout << sol.serverCluster[i] << std::endl;*/
}

double getElapsedSeconds(const std::chrono::steady_clock::time_point& start_time)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
}
