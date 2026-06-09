#pragma once
#ifndef _COMMON_FUNC_H
#define _COMMON_FUNC_H
#include"type_def.h"

int computeMoveDelta(const KMBCPInstance& G, const Solution& sol, int si, int fromK, int toK);
int computeMoveDelta_improved(const KMBCPInstance& G, const Solution& sol, const Help_struc_move& struc_move, int si, int fromK, int toK);
int computeSwapDelta(const KMBCPInstance& G, const Solution& sol, int si, int sj);
int computeSwapDelta_improved(const KMBCPInstance& G, const Solution& sol, const Help_struc_move& struc_move, int si, int sj);
void applyMove(const KMBCPInstance& G, Solution& sol, Help_struc_move& struc_move, int si, int toK, int best_delta);
void applySwap(const KMBCPInstance& G, Solution& sol, Help_struc_move& struc_move, int si, int sj, int best_delta);

void generateInitialSolutionRandom(const KMBCPInstance& G, Solution& sol, Help_struc_move& struc_move, int seed);
void copy_solution(const Solution& sol_src, Solution& sol_des);
void check_solution(const Solution& sol, const KMBCPInstance& G, int cost);
void output_res(const std::string outFileName, const std::string input_file, const Solution& best_sol, double first_found_time,
	double miningTime,
	int seed);

Solution tabuSearch(const KMBCPInstance& ins, Solution& current_sol, int maxIter, int tabuTenure, std::mt19937& gen);
Solution solution_based_tabuSearch(const KMBCPInstance& ins, Solution& current_sol, int maxIter, std::mt19937& gen, std::mt19937_64& gen_64);

void pattern_mining(const KMBCPInstance& instance, std::vector<Pattern>& mp_res);

void evolutionary_search(const KMBCPInstance& ins, const std::vector<Pattern>& mp_res, std::mt19937& gen, std::mt19937_64& gen_64);

/*elapsed time in seconds*/
double getElapsedSeconds(const std::chrono::steady_clock::time_point& start_time);

#endif