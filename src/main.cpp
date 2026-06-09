#include"type_def.h"
#include"common_func.h"

Variable_gb Vars_gb;
Paras_alg Paras_gb;

int main(int argc, char* argv[]) 
{
    try {
        
        if (argc < 4) {
            std::cout << "Alg usage: input_file random_seed output_sol";
            exit(-1);
        }
        //std::string input_file = "F:\\clique problem\\K-CMBCP(k-clustering minimum biclique completion problem)\\Data_KCluster\\DataK-clustering\\150-150-10-0.3.dat";
        std::string input_file = argv[1];
        int seed = atoi(argv[2]);      
        //std::string output_sol = "E:\\output\\PMHA_sol.txt";
        std::string output_sol = argv[3];

        KMBCPInstance ins;
        std::vector<Pattern> mp_res;             //maxima patterns results
        std::mt19937 gen(seed);
        std::mt19937_64 gen_64(seed);

        Vars_gb.start_time = std::chrono::steady_clock::now();          //以后不能用clock函数计算时间了，太不准确

        ins.readKMBCPInstance(input_file);
        //printKMBCPInstance("F:\\clique problem\\K-CMBCP(k-clustering minimum biclique completion problem)\\Data_KCluster\\DataK-clustering\\15-15-5-0.5.dat", ins);
        
        int p_size = 50;                         //population size
        double mu_ratio = 0.10;                  //mutatant ratio        
        int maxIter_ts = 10000;                  //maximum consecutive iterations of non-improvement for tabu search
        int tabuTenure = 15;                     //tabu tenure for the attibute-based tabu search        
        int min_sup = static_cast<int>(ins.n * 0.4);    //min support for the pattern mining
        Paras_gb.set_Alg_paras(p_size, mu_ratio, maxIter_ts, tabuTenure, min_sup);

        if (ins.n <= 100)
            Vars_gb.time_limit = 300;
        else
            Vars_gb.time_limit = 1800;

        pattern_mining(ins, mp_res);
        double miningTime = getElapsedSeconds(Vars_gb.start_time);
        std::cout << "mp_res.size=" << mp_res.size() << std::endl;

        evolutionary_search(ins, mp_res, gen, gen_64);
        output_res(output_sol, input_file, Vars_gb.Best_sol, Vars_gb.first_time, miningTime, seed);

        auto now = std::chrono::steady_clock::now();
        double elapsed_time = std::chrono::duration<double>(now - Vars_gb.start_time).count();
        std::cout << std::endl << "the evolutionary search is finisehd, input_file = " << input_file << ", best_v = " << Vars_gb.Best_sol.cost << ", best_time = " <<
            Vars_gb.first_time << ", random_seed = " << seed << ", elapsed_time = " << elapsed_time << std::endl;               
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}

