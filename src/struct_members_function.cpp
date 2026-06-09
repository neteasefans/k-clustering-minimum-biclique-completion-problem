#include"type_def.h"
#include <fstream>


// 读取K-MBCP算例
void KMBCPInstance::readKMBCPInstance(const std::string & filename)
{
    std::ifstream infile(filename);
    std::string line;

    if (!infile.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    // 读取第一行：n m K
    if (!std::getline(infile, line)) {
        throw std::runtime_error("Empty file or cannot read first line: " + filename);
    }
    if (!line.empty() && line.back() == '\r')
        line.pop_back();                                 // 去掉 Windows 回车

    std::cout << "filename=" << filename << std::endl;
    std::istringstream iss(line);
    if (!(iss >> this->n >> this->m >> this->K)) {
        throw std::runtime_error("Invalid first line format: " + line);
    }

    // 初始化 edges
    this->edges.assign(static_cast<size_t>(this->n), std::vector<bool>(static_cast<size_t>(this->n) + static_cast<size_t>(this->m), false));

    // 读取边数据
    while (std::getline(infile, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();                         // 去掉回车

        if (line.empty())
            continue;                               // 跳过空行

        std::stringstream ss(line);
        size_t u, v;
        if (!(ss >> u >> v)) {
            std::cerr << "Illegal edge line, skipping: " << line << std::endl;
            continue;
        }
        if (u < 1 || u > this->n || v < 1 || v > static_cast<size_t>(this->n) + static_cast<size_t>(this->m)) {
            std::cerr << "Edge out of range, skipping: " << line << std::endl;
            continue;
        }
        this->edges[u - 1][v - 1] = true;
    }

    // 构建 serverToClients 和 clientToServers
    this->serverToClients.assign(static_cast<size_t>(this->n), std::vector<int>());
    this->clientToServers.assign(static_cast<size_t>(this->n) + static_cast<size_t>(this->m), std::vector<int>());

    for (int i = 0; i < this->n; ++i) {
        for (int j = this->n; j < this->n + this->m; ++j) {
            if (this->edges[i][j]) {
                this->serverToClients[i].push_back(j);
                this->clientToServers[j].push_back(i);
            }
        }
    }
    infile.close();
}

/*打印K-MBCP算例*/
void KMBCPInstance::printKMBCPInstance(const std::string& filename)
{
    std::cout << "Instance read: n=" << this->n << ", m=" << this->m << ", K=" << this->K << std::endl;
    for (int i = 0; i < this->n; i++) {
        for (size_t j = 0; j < this->edges[i].size(); j++) {
            if (this->edges[i][j])
                std::cout << i + 1 << " " << j + 1 << std::endl;
        }
    }
    std::cout << "segment line -----------------------" << std::endl;
    std::cout << "The following is adjacent matrix----" << std::endl;
    for (int i = 0; i < this->n; i++) {
        std::cout << "i=" << i << std::endl;
        for (size_t j = 0; j < this->serverToClients[i].size(); j++)
            std::cout << this->serverToClients[i][j] << std::endl;
    }
    std::cout << "print instance finisehd" << std::endl;
}


/*calculate the cost of a solution
* time complexity: O(m * K)
*/
void Solution::set_cost(const KMBCPInstance& G)
{
    this->cost = 0;
    for (int j = G.n; j < G.n + G.m; j++) {
        for (int k = 0; k < G.K; k++) {
            if (this->clientDegreeInCluster[j][k] > 0)
                this->cost += this->sizeCluster[k] - this->clientDegreeInCluster[j][k];
        }
    }
}


/*构建solution的成员， 根据cluster assignment*/
void Solution::build_solution_data(const KMBCPInstance& G, const std::vector<int>& clusters_of)
{
    this->reset(G.n, G.m, G.K);
    this->serverCluster = clusters_of;
    for (int i = 0; i < G.n; i++) {
        int k = this->serverCluster[i];
        for (int j : G.serverToClients[i]) {
            this->clientDegreeInCluster[j][k]++;
        }
        this->servers_in_cluster[k].push_back(i);
    }
    for (int j = G.n; j < G.n + G.m; j++) {
        for (int k = 0; k < G.K; k++) {
            if (this->clientDegreeInCluster[j][k] > 0)
                this->clients_in_cluster[k].push_back(j);
        }
    }
    for (int i = 0; i < G.n; i++) {
        if (this->serverCluster[i] >= 0 && this->serverCluster[i] < G.K)
            this->sizeCluster[this->serverCluster[i]]++;
    }
    set_cost(G);
}

/*初始化move helping structure, for fast incremental evaluation
* complexity: O(K * m + n*k*d(i))
*/
void Help_struc_move::build_strucmove_data(const KMBCPInstance& ins, const Solution& sol)
{

    for (int k = 0; k < ins.K; k++) {
        for (int j : sol.clients_in_cluster[k]) {
            if (sol.clientDegreeInCluster[j][k] == 1)
                this->servers_eq_one_clients[k].push_back(j);
            else if (sol.clientDegreeInCluster[j][k] > 1)
                this->servers_ge_one_clients[k].push_back(j);
        }
    }

    for (int i = 0; i < ins.n; i++) {
        for (int j : ins.serverToClients[i]) {
            bool exists = false;
            for (int k = 0; k < ins.K; k++) {
                exists = false;
                if (sol.clientDegreeInCluster[j][k] == 1) {
                    this->servers_eq_one_intersect_clients[i][k].push_back(j);
                }
                else if (sol.clientDegreeInCluster[j][k] > 1)
                    this->servers_ge_one_intersect_clients[i][k].push_back(j);
                for (int j2 : sol.clients_in_cluster[k]) {
                    if (j == j2) {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    this->server_adjs_not_in_clients[i][k].push_back(j);
            }
        }
    }
    for (int i = 0; i < ins.n; i++) {
        for (int k = 0; k < ins.K; k++) {
            for (int j : sol.clients_in_cluster[k]) {
                if (ins.edges[i][j])
                    this->server_adjs_in_clients[i][k].push_back(j);
            }
        }
    }

}

/*set parameters of the algorithm*/
void Paras_alg::set_Alg_paras(int p_size, double mu_ratio, int max_iter, int tabu_tenure, int min_sup)
{
    this->pop_size = p_size;
    this->mutate_ratio = mu_ratio;
    this->maxIter_ts = max_iter;
    this->tabuTenure = tabu_tenure;
    this->min_sup = min_sup;
}