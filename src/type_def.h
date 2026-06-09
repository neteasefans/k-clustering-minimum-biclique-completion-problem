#pragma once
#ifndef _TYPE_DEF_H
#define _TYPE_DEF_H
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm> 
#include <random>
#include<climits>
#include<chrono>

/*算例结构*/
struct KMBCPInstance {
    int n; // 左部点数量
    int m; // 右部点数量
    int K; // 簇的数量
    std::vector<std::vector<bool>> edges;                // n x (n+m) 的边集， 表示两个顶点是否有边相连   
    std::vector<std::vector<int>> serverToClients;       //邻接矩阵， 记录每个左部点（服务点）有边相连的右边点（客户点）集合
    std::vector<std::vector<int>> clientToServers;       //邻接矩阵， 记录每个右部点（客户点）有边相连的左边点（服务点）集合

    void readKMBCPInstance(const std::string& filename);               //read KMBCP instance
    void printKMBCPInstance(const std::string& filename);        //print KMBCP instance
};

/*解结构*/
struct Solution {
    std::vector<int> serverCluster;                       // serverCluster[i] = k 表示服务点i属于第k个簇
    std::vector<std::vector<int>> clientDegreeInCluster;  // clientDegreeInCluster[j][k] 表示客户端j在第k个簇内的邻接服务点个数， 若==0， 则表示客户点j不在该cluster中，这个概念非常重要， chatgpt给的提示太厉害
    std::vector<int> sizeCluster;                         //sizeCluster[k]表示cluster k中服务点个数
    std::vector<std::vector<int>> servers_in_cluster;     //servers_in_cluster[k] 表示cluster k中服务点集合
    std::vector<std::vector<int>> clients_in_cluster;     //clients_in_cluster[k] 表示cluster k中客户点集合
    long long cost;                                         // 当前解的目标值

     // 默认构造函数：空
    Solution() : cost(0) {}

    // 参数构造函数：分配空间并初始化为0
    Solution(int n, int m, int K) : cost(0) {
        reset(n, m, K);
    }

    // 重置函数：清空并重新分配
    void reset(int n, int m, int K) {
        cost = 0;
        serverCluster.assign(n, -1);                            // 初始为 -1 表示未分配
        clientDegreeInCluster.assign(size_t(n) + size_t(m), std::vector<int>(K, 0));
        sizeCluster.assign(K, 0);
        servers_in_cluster.assign(K, {});
        clients_in_cluster.assign(K, {});
    } 

    void build_solution_data(const KMBCPInstance& G, const std::vector<int>& clusters_of);      //initializing the solutin data, according to the input cluster assignment

    void set_cost(const KMBCPInstance& G);                                                  //calculate the solution cost
};

/*简化的解结构*/
struct Individual {
    std::vector<int> cluster;                                // serverCluster[i] = k 表示服务点i属于第k个簇
    int fitness = 0;                                         // 当前解的目标值

    //默认构造函数
    Individual(int n)                                        // 初始为 -1 表示未分配
        : cluster(n, -1)
    {}
};


/*for fast incremental evaluation
* for one-move and swap move operators
*/
struct Help_struc_move
{    
    std::vector<std::vector<int>> servers_ge_one_clients;     //servers_ge_one_clients[k]表示cluster k中，deg > 1的客户点集合
    std::vector<std::vector<int>> servers_eq_one_clients;     //servers_eq_one_clients[k]表示cluster k中，deg = 1的客户点集合
    std::vector<std::vector<std::vector<int>>> servers_eq_one_intersect_clients;     //servers_eq_one_intersect_clients[i][k]表示cluster k中，deg = 1, 且与服务点i相邻的客户点集合
    std::vector<std::vector<std::vector<int>>> servers_ge_one_intersect_clients;     //servers_ge_one_intersect_clients[i][k]表示cluster k中，deg > 1, 且与服务点i相邻的客户点集合

    //std::vector<std::vector<std::vector<std::vector<int>>>> servers_eq_one_non_intersect2_clients;  //servers_eq_one_non_intersect2_clients[si][sj][k]表示cluster k中，deg = 1, 且与服务点si、sj均不相邻的客户点集合
    //std::vector<std::vector<std::vector<std::vector<int>>>> servers_eq_one_intersect2_clients;     //servers_eq_one_intersect2_clients[si][sj][k]表示cluster k中，deg = 1, 且与服务点si、sj都相邻的客户点集合
    
    std::vector<std::vector<std::vector<int>>> server_adjs_not_in_clients;     //server_adjs_not_in_clients[i][k]表示与服务点i相邻，且不在cluster k中的客户点集合
    std::vector<std::vector<std::vector<int>>> server_adjs_in_clients;          //server_adjs_in_clients[i][k]表示与服务点i相邻，且在cluster k中的客户点集合

    //参数构造函数
    Help_struc_move(int n, int K) {
        servers_eq_one_clients.assign(K, std::vector<int>());
        servers_ge_one_clients.assign(K, std::vector<int>());
        servers_eq_one_intersect_clients.assign(n, std::vector<std::vector<int>>(K));
        servers_ge_one_intersect_clients.assign(n, std::vector<std::vector<int>>(K));
        server_adjs_not_in_clients.assign(n, std::vector<std::vector<int>>(K));
        server_adjs_in_clients.assign(n, std::vector<std::vector<int>>(K));
    }

    void build_strucmove_data(const KMBCPInstance& ins, const Solution& sol);      //initializing the strucmove_data according the input solution
};

/*
* for pattern mining
*/
struct Pattern
{
    std::vector<int> items;
    std::vector<int> tidset;
    int support;
};

/*
* globla variables
*/
struct Variable_gb
{
    Solution Best_sol;                          //算法找到的最好解
    double first_time;                          //找到最好解，所需的时间
    std::chrono::time_point<std::chrono::steady_clock> start_time;                          //the starting time for the algorithm
    double time_limit;                          //time limit for the algorithm
};

/*
* parameters for the algorithm
*/
struct Paras_alg
{
    //parameters for evolutionary search
    int pop_size;                //population size
    double mutate_ratio;         //mutation ratio
    double weight_alpha;         //unused, the weight of fitnees and , diversity
    int generations;             //unused,  the number of generations

    //parameters for tabu search
    int maxIter_ts;             //maximum number of non-improvement  consecutive iterations
    int tabuTenure;             //tabu tenure
  
    //parameters for pattern mining
    int min_sup;                                //minimal support for pattern mining

    void set_Alg_paras(int p_size, double mu_ratio, int max_iter, int tabu_tenure, int min_sup);                   //set the algorithm parameters
};

extern Variable_gb Vars_gb;
extern Paras_alg Paras_gb;
#endif