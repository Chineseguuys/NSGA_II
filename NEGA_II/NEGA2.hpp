#pragma once

#ifndef __NEGA2__
#define __NEGA2__

#ifndef _IOSTREAM_
#include <iostream>
#include <fstream>
#endif // !_IOSTREAM_

#ifndef _VECTOR_
#include <vector>
#endif // !_VECTOR_

#ifndef _STRING_
#include <string>
#endif // !_STRING_

#ifndef _CTIME_
#include <ctime>
#endif // !_CTIME_

#ifndef _UTILITY_
#include <utility>
#endif // !_UTILITY_

#ifndef _ALGORITHM_
#include<algorithm>
#endif // !_ALGORITHM_


#ifndef _INF_
#define INF  1e+14     //位于边缘的解向量的距离值为正无穷
#endif // !_INF_

#ifndef _EPS_
#define EPS  1e-14     //定义一个无穷小量， 方便进行double类型的大小比较
#endif // !_EPS_


using namespace std;

class Random_Integer {
	public:
		int get_RandomInteger(int lower_limit, int upper_limit)              //返回一定范围内的整数
		{
			return (rand() % (upper_limit - lower_limit + 1) + lower_limit);
		}

		double get_randomReal(double lower_limits, double  upper_limits)      //返回一定范围的随机实数
		{
			return (lower_limits + (upper_limits - lower_limits)*rand() / (RAND_MAX + 1));
		} 

		Random_Integer()
		{
			srand(unsigned int(time(NULL)));
		}

		Random_Integer(unsigned int seed)
		{
			srand(seed);
		}
};


typedef void (*functionType)(double *, double *);             // 所有要测试的函数的函数指针


class Individual         //每一个解向量类
{
	private:
		int Number_Variables;    //变量的个数
		int Number_ObjectFunction;   //测试函数当中相互约束的函数的个数(需要综合考虑的函数的个数)
		int be_dominated_number;     //被多少个个体所支配，将会用于排名(Rank)，初始化为零
		int Serial_Number;
		int Rank;    //我的排名等级是多少
		bool Evaluated;
		static std::vector< std::pair<double, double> > Value_Limits;    //每一个变量的作用的范围
		static Random_Integer random_generator;

	public:
		/*
			经常访问的变量
		*/
		double Crowded_Distance;  //拥挤距离
		double *Variables;
		double *ObjectFunction_Value;   //每一个约束函数的值
		bool Has_Ranked;
		bool Has_Crowded;
		std::vector<Individual*> Dominating_Inividual;                   // 支配了哪一些变量,以指针的形式指向其支配的对象
		static int Get_Obj;

		Individual(int number_variables, 
				   int number_objectfuncton,
				   int serial_number
				   )
		{
			this->Number_Variables = number_variables;
			this->Number_ObjectFunction = number_objectfuncton;
			this->be_dominated_number = 0;
			this->Serial_Number = serial_number;
			this->Rank = 0;
			this->Crowded_Distance = 0.0;
			this->Evaluated = false;
			this->Variables = new double[this->Number_Variables];
			this->Has_Ranked = false;
			this->Has_Crowded = false;
			this->ObjectFunction_Value = new double[this->Number_ObjectFunction];
		}

		~Individual()
		{
			delete this->Variables;
			delete this->ObjectFunction_Value;
			this->Variables = NULL;
			this->ObjectFunction_Value = NULL;
			//printf_s("delete indivdual\n");
		}

		void operator=(const Individual& a)
		{
			this->Crowded_Distance = a.Crowded_Distance;
			this->Number_Variables = a.Number_Variables;
			this->Number_ObjectFunction = a.Number_ObjectFunction;
			this->Rank = a.Rank;
			this->be_dominated_number = a.be_dominated_number;
			this->Serial_Number = a.Serial_Number;
			this->Evaluated = a.Evaluated;
			for (int i = 0; i < this->Number_Variables; i++)
				this->Variables[i] = a.Variables[i];
			for (int i = 0; i < this->Number_ObjectFunction; i++)
				this->ObjectFunction_Value[i] = a.ObjectFunction_Value[i];
		}

		Individual(const Individual& a)
		{
			this->Crowded_Distance = a.Crowded_Distance;
			this->Number_Variables = a.Number_Variables;
			this->Number_ObjectFunction = a.Number_ObjectFunction;
			this->Rank = a.Rank;
			this->be_dominated_number = a.be_dominated_number;
			this->Serial_Number = a.Serial_Number;
			this->Evaluated = a.Evaluated;
			this->Variables = new double[this->Number_Variables];
			this->ObjectFunction_Value = new double[this->Number_ObjectFunction];
			for (int i = 0; i < this->Number_Variables; i++)
				this->Variables[i] = a.Variables[i];
			for (int i = 0; i < this->Number_ObjectFunction; i++)
				this->ObjectFunction_Value[i] = a.ObjectFunction_Value[i];
		}

		void init_solutionVector()           //初始化每一个解向量
		{
			double lower_limit;
			double upper_limit;
			for (int i = 0; i < this->Number_Variables; i++)
			{
				lower_limit = Value_Limits[i].first;
				upper_limit = Value_Limits[i].second;
				this->Variables[i] = Individual::random_generator.get_randomReal(lower_limit, upper_limit);
			}
		}

		void calculate_function_value(functionType function)
		{
			function(this->Variables, this->ObjectFunction_Value);
			this->Evaluated = true;
		}

		double* get_value()
		{
			return this->ObjectFunction_Value;
		}

		void append_dominating(Individual* dominated)         //指向被其统治的个体，并增加统治数量
		{
			this->Dominating_Inividual.push_back(dominated);
		}

		void reset_evaluated()          //重置其计算状态，在变异之后可能需要重新计算目标函数的值
		{
			this->Evaluated = false;
		}

		bool get_evaluatedState()
		{
			return this->Evaluated;
		}

		void change_be_dominatedNumber(int num)
		{
			this->be_dominated_number = this->be_dominated_number + num;
		}

		int get_dominatedNum()
		{
			return this->be_dominated_number;
		}

		void reset_dominatedNum()
		{
			this->be_dominated_number = 0;
			this->Dominating_Inividual.clear();
			this->Has_Ranked = false;
		}

		void set_rank(int rank)
		{
			this->Rank = rank;
			this->Has_Ranked = true;
		}

		void reduce_domulated()
		{
			this->be_dominated_number = this->be_dominated_number - 1;
		}

		int get_rank()
		{
			return this->Rank;
		}

		int get_serial_num()
		{
			return this->Serial_Number;
		}

		static void set_valuelimits(const std::vector< std::pair<double, double> >& value_limits)
		{
			Individual::Value_Limits.assign(value_limits.begin(), value_limits.end());
		}
};

// static 变量的声明
std::vector< std::pair<double, double> > Individual::Value_Limits;    //占用的空间比较的大，设置为 static 变量比较好
Random_Integer Individual::random_generator;
int Individual::Get_Obj = 0;


// compare 函数
bool comp(Individual* a,Individual* b)
{
	return (a->get_rank() < b->get_rank());
}

bool comp_obj(Individual* a, Individual* b)        //按照函数值从大到小排列
{
	return (a->ObjectFunction_Value[Individual::Get_Obj] > b->ObjectFunction_Value[Individual::Get_Obj]);
}

bool comp_dist(Individual* a, Individual* b)        //按照函数值从小到大排列
{
	return (a->Crowded_Distance > b->Crowded_Distance);
}

//NSGA2 类
class NSGA2
{
	private:
		int Number_Variables;
		int Number_ObjectFunction;
		int Population_Size;
		double Mutation_Rate;                  // 变异的概率
		double Crossover_Rate;                 //交叉的概率
		double Cross_DistrIndex;             //交叉分布指数
		double Mutate_DistrIndex;           //变异的交叉分布指数
		int Rank_Num;
		std::vector<Individual*> population;                      //每一代中的种群
		std::vector<Individual*> prerant_chirdern;        // 父代和子代的综合种群
		functionType testing_function;                // 用于测试的函数
		std::vector< std::pair<double, double> > limits_realvar;        //每一个变量的作用范围
		static Random_Integer NSGA2_random_generator;                   //是不是static都不所谓，反正只产生一个实例

	public:
		NSGA2(int number_variables, int number_objectfunction, int population_size, 
			  double mutation_rate, double crossover_rate, double cross_distrindex, double mutate_distrindex,
			  functionType function)
		{
			this->Number_Variables = number_variables;
			this->Number_ObjectFunction = number_objectfunction;
			this->Population_Size = population_size;
			this->Mutation_Rate = mutation_rate;
			this->Crossover_Rate = crossover_rate;
			this->Cross_DistrIndex = cross_distrindex;
			this->Mutate_DistrIndex = mutate_distrindex;
			this->Rank_Num = 0;
			this->testing_function = function;
		}

		~NSGA2()
		{
			for (std::vector<Individual*>::iterator it = population.begin(); it != population.end(); it++)
			{
				if (NULL != *it)
				{
					delete *it;
					*it = NULL;
				}
			}
		}

		void set_limits(double value_limits[][2], int size)           //设置每一个变量的取值范围
		{
			if (size != this->Number_Variables)
			{
				std::cout << "You have gave a wrong size of value limits" << '\n';
				return;
			}
			std::pair<double, double> temp;
			for (int i = 0; i < size; i++)
			{
				temp.first = value_limits[i][0];
				temp.second = value_limits[i][1];
				limits_realvar.push_back(temp);
			}
		}

		void initial_population();               //初始化种群中的每一个解向量
		void deciding_domination(Individual* a, Individual* b);     // 确定种群中的两个个体的支配关系 (a支配b,或b支配a,或不相互支配)
		void simulated_binary_x(const Individual& a, const Individual& b, Individual& child_a, Individual& child_b);       //使用模拟二进制交叉，和多项式变异来产生两个子代的个体
		void polynomial_mutation(Individual& a);                       //多项式变异
		void nondominated_sort(std::vector<Individual*>& pop);          //非支配排序
		void iteration_rank(Individual* individual, int rank);         // rank 
		void calculate_croweddist(std::vector<Individual*>& pop);         //计算拥挤距离
		void evoluate(int iteration_times);
		void run(int iteration_times);
		double get_CDI();                                              // 获取交叉分布因子
		void display();
}; 
// static 变量的声明
Random_Integer NSGA2::NSGA2_random_generator;


void NSGA2::initial_population()         //初始化了所有的解向量，并计算了所有的目标函数的值
{
	for (int i = 0; i < this->Population_Size; i++)
	{
		Individual* temp = new Individual(this->Number_Variables, this->Number_ObjectFunction, i + 1);
		this->population.push_back(temp);
	}
	Individual::set_valuelimits(this->limits_realvar);
	for (int i = 0; i < this->Population_Size; i++)
	{
		this->population[i]->init_solutionVector();          //初始化所有的解向量的值
		this->population[i]->calculate_function_value(this->testing_function);  //计算所有的目标函数的值
	}
	this->nondominated_sort(this->population);
	this->calculate_croweddist(this->population);
}

void NSGA2::deciding_domination(Individual* a, Individual* b)    //确定两个变量之间的支配关系
{
	/*
		支配关系一定是支配或者是非支配，不需要考虑完全相等的情况
	*/
	//先判断是否计算了函数值
	a->calculate_function_value(this->testing_function);
	b->calculate_function_value(this->testing_function);

	double *value_a = a->get_value();
	double *value_b = b->get_value();
	int count = 0;
	for (int i = 0; i < this->Number_ObjectFunction; i++)
	{
		if (value_a[i] < value_b[i])
			count += 1;
		else
			break;
	}
	if (count == this->Number_ObjectFunction)
	{
		a->append_dominating(b);
		b->change_be_dominatedNumber(1);
	}
	count = 0;
	for (int i = 0; i < this->Number_ObjectFunction; i++)
	{
		if (value_b[i] < value_a[i])
			count += 1;
		else
			break;
	}
	if (count == this->Number_ObjectFunction)
	{
		b->append_dominating(a);
		a->change_be_dominatedNumber(1);
	}
	return;
}

double NSGA2::get_CDI()
{
	/*
		得到NSGA2中的交叉分布指数
	*/
	return this->Cross_DistrIndex;
}

void NSGA2::simulated_binary_x(const Individual& a, const Individual& b, Individual& child_a, Individual& child_b)
{
	/*
		模拟二进制交叉算法 (SBX)
		注意：这里产生的子代并没有计算具体的函数值，默认的函数值都是零(或者任意数)
	*/
	double v_low, v_up, v_low_limit, v_up_limit;
	double  u;         //0 到 1 之间的随机数
	double alpha_1, alpha_2;
	double bata_1, bata_2;
	double xo1, xo2;

	if (NSGA2::NSGA2_random_generator.get_randomReal(0.0, 1.0) < this->Crossover_Rate)
	{
		for (int i = 0; i < this->Number_Variables; i++)
		{
			if (fabs(a.Variables[i] - b.Variables[i] > EPS))
			{
				if (a.Variables[i] < b.Variables[i])     //确定大小关系，方便后面的数学计算
				{
					v_low = a.Variables[i];
					v_up = b.Variables[i];
				}
				else
				{
					v_low = b.Variables[i];
					v_up = a.Variables[i];
				}
				v_low_limit = this->limits_realvar[i].first;     //确定上下限
				v_up_limit = this->limits_realvar[i].second;
				alpha_1 = pow((1 + 2 * (v_low - v_low_limit) / (v_up - v_low)), -(this->Cross_DistrIndex + 1));
				alpha_1 = 2 - alpha_1;
				alpha_2 = pow((1 + 2 * (v_up_limit - v_up) / (v_up - v_low)), -(this->Cross_DistrIndex + 1));
				alpha_2 = 2 - alpha_2;
				u = NSGA2::NSGA2_random_generator.get_randomReal(0.0, 1.0);

				if (u <= (1 / alpha_1))
				{
					bata_1 = pow((u*alpha_1), (1 / (this->Cross_DistrIndex + 1)));
				}
				else
				{
					bata_1 = 1 / (2 - u*alpha_1);
					bata_1 = pow(bata_1, (1 / (this->Cross_DistrIndex + 1)));
				}
				if (u <= (1 / alpha_2))
				{
					bata_2 = u*alpha_2;
					bata_2 = pow(bata_2, (1 / (this->Cross_DistrIndex + 1)));
				}
				else
				{
					bata_2 = 1 / (2 - u*alpha_2);
					bata_2 = pow(bata_2, (1 / (this->Cross_DistrIndex + 1)));
				}

				xo1 = 0.5*(v_low + v_up - bata_1*(v_up - v_low));
				xo2 = 0.5*(v_low + v_up + bata_2*(v_up - v_low));

				// 值给哪个儿子是否很重要
				xo1 = min(max(xo1, v_low_limit), v_up_limit);
				xo2 = min(max(xo2, v_low_limit), v_up_limit);
				//printf_s("SBX :(%f, %f)\n", xo1, xo2);
				child_a.Variables[i] = xo1;
				child_b.Variables[i] = xo2;
			}
			else
			{
				/*
					某一位的值几乎相同，直接传递给下一代
				*/
				child_a.Variables[i] = a.Variables[i];
				child_b.Variables[i] = b.Variables[i];
				//printf_s("SBX :(%f, %f)\n", a.Variables[i], b.Variables[i]);
			}
		}
	}
	else
	{
		/*
		     父代直接继承给子代个体， 当不进行交叉的时候
		*/
		child_a = a;
		child_b = b;
	}
	child_a.reset_evaluated();
	child_b.reset_evaluated();
}

void NSGA2::polynomial_mutation(Individual& a)            //多项式变异
{
	double theta_q;
	double x_p;
	double x_o;
	double u;      //0 到 1 之间的随机变量
	double v_low_limit, v_up_limit;

	for (int i = 0; i < this->Number_Variables; i++)
	{
		if (NSGA2::NSGA2_random_generator.get_randomReal(0.0, 1.0) < this->Mutation_Rate)
		{
			x_p = a.Variables[i];
			v_low_limit = this->limits_realvar[i].first;
			v_up_limit = this->limits_realvar[i].second;
			u = NSGA2::NSGA2_random_generator.get_randomReal(0.0, 1.0);
			if (u <= 0.5)
			{
				double xy = 1.0 - (x_p - v_low_limit) / (v_up_limit - v_low_limit);
				theta_q = 2.0*u + (1 - 2.0*u)*pow(xy, (this->Mutate_DistrIndex + 1));
				theta_q = pow(theta_q, (1.0 / (this->Mutate_DistrIndex + 1))) - 1.0;
			}
			else
			{
				/*theta_q = 2 * (1 - u) + 2 * (u - 0.5)*(1 - (v_up_limit - x_p) / (v_up_limit - v_low_limit));
				theta_q = pow(theta_q, (1 / (this->Mutate_DistrIndex + 1)));
				theta_q = 1 - theta_q;*/
				double xy = 1 - (v_up_limit - x_p) / (v_up_limit - v_low_limit);
				theta_q = 2.0*(1.0 - u) + 2.0*(u - 0.5)*pow(xy, (this->Mutate_DistrIndex + 1));
				theta_q = 1.0 - pow(theta_q, (1 / (this->Mutate_DistrIndex + 1)));
			}
			x_o = x_p + theta_q*(v_up_limit - v_low_limit);
			//cout << "mutation : " << a.Variables[i] << ", " << x_o << '\n';
			//需不需要判断 x_o 的越界的情况    不太清楚
			if (x_o < v_low_limit)
				x_o = v_low_limit;
			if (x_o > v_up_limit)
				x_o = v_up_limit;
			
			a.Variables[i] = x_o;
		}
		else
		{
			continue;
		}
	}
	a.reset_evaluated();
}

void NSGA2::nondominated_sort(std::vector<Individual*>& pop)
{
	//先计算整个种群中的非支配关系
	for (std::vector<Individual*>::iterator it = pop.begin(); it != pop.end(); it++)
	{
		(*it)->reset_dominatedNum();
	}
	for (int i = 0; i < pop.size() - 1; i++)
	{
		for (int j = i+1; j < pop.size(); j++)
		{
			this->deciding_domination(pop[i], pop[j]);
		}
	}
	//重写快速非支配排序
	for (std::vector<Individual*>::iterator it = pop.begin(); it != pop.end(); it++)
	{
		if ((*it)->get_dominatedNum() == 0)
			this->iteration_rank(*it, 1);
	}

	sort(pop.begin(), pop.end(), comp);   //根据rank排序
}


void NSGA2::iteration_rank(Individual* a, int rank)    //循环迭代的过程
{
	if ((a->get_dominatedNum() == 0) && (a->Has_Ranked == false))
	{
		a->set_rank(rank);
		if (a->Dominating_Inividual.size() > 0)
		{
			int size = a->Dominating_Inividual.size();
			for (int i = 0; i < size; i++)
			{
				a->Dominating_Inividual[i]->reduce_domulated();
				this->iteration_rank(a->Dominating_Inividual[i], rank + 1);
			}
		}
	}
	else
	{
		return;
	}
}

void NSGA2::calculate_croweddist(std::vector<Individual*>& pop)    // 计算拥挤距离，并根据拥挤距离进行排序
{
	//cout << "calculate_croweddist" << '\n';
	for (std::vector<Individual*>::iterator it = pop.begin(); it != pop.end(); it++)
	{
		(*it)->Crowded_Distance = 0.0;
		(*it)->Has_Crowded = false;
	}
	std::vector<Individual*>::iterator it = pop.end();
	it--;
	this->Rank_Num = (*it)->get_rank();
	std::vector<Individual*>::iterator start = pop.begin();
	std::vector<Individual*>::iterator end = pop.begin();
	int count = 0;
	double fmax, fmin;
	for (int i = 1; i <= this->Rank_Num; i++)
	{
		count = 0;
		while ((*end)->get_rank() == i)
		{
			count++;
			++end;
			if (end == pop.end())
			{
				break;
			}
		}
		if (count == 1)
		{
			(*start)->Crowded_Distance = INF;
			(*start)->Has_Crowded = true;
		}
		else if (count == 2)
		{
			(*start)->Crowded_Distance = INF;
			(*start)->Has_Crowded = true;
			(*(end - 1))->Crowded_Distance = INF;
			(*(end - 1))->Has_Crowded = true;
		}
		else
		{
			for (int j = 0; j < this->Number_ObjectFunction; j++)
			{
				Individual::Get_Obj = j;
					sort(start, end, comp_obj);
					fmax = (*start)->ObjectFunction_Value[j];
					fmin = (*(end - 1))->ObjectFunction_Value[j];
					(*start)->Crowded_Distance = INF;
					(*start)->Has_Crowded = true;
					(*(end - 1))->Crowded_Distance = INF;
					(*(end - 1))->Has_Crowded = true;
					for (int k = 1; k < count - 1; k++)
					{
						if ((*(start + k))->Crowded_Distance == INF)
							continue;
						else
						{
							//cout << "fabs : "<<j <<"-" << fabs((*(start + k - 1))->ObjectFunction_Value[j] - (*(start + k + 1))->ObjectFunction_Value[j]) << '\n';
							(*(start + k))->Crowded_Distance += (fabs((*(start + k - 1))->ObjectFunction_Value[j] - (*(start + k + 1))->ObjectFunction_Value[j]) / (fmax - fmin));
							(*(start + k))->Has_Crowded = true;
						}
					}
			}
		}
		sort(start, end, comp_dist);       //根据拥挤距离排序
		start = end;
	}
}



void NSGA2::evoluate(int iteration_times)
{
	std::vector<Individual*> temp;
	//二元锦标赛
	for (int i = 0; i < this->Population_Size; i++)
	{
		int a = NSGA2::NSGA2_random_generator.get_RandomInteger(0, this->Population_Size - 1);
		int b = NSGA2::NSGA2_random_generator.get_RandomInteger(0, this->Population_Size - 1);
		if (this->population[a]->get_rank() < this->population[b]->get_rank())
		{
			Individual* one = new Individual(*(this->population[a]));
			temp.push_back(one);
		}
		else if (this->population[b]->get_rank() < this->population[a]->get_rank())
		{
			Individual* one = new Individual(*(this->population[b]));
			temp.push_back(one);
		}
		else
		{
			if (this->population[a]->Crowded_Distance < this->population[b]->Crowded_Distance)
			{
				Individual* one = new Individual(*(this->population[b]));
				temp.push_back(one);
			}
			else
			{
				Individual* one = new Individual(*(this->population[a]));
				temp.push_back(one);
			}
		}
	}
	std::swap(this->population, temp);
	this->prerant_chirdern.clear();
	//交叉变异
	for (int i = 0; int(i < this->Population_Size / 2); i++)
	{
		int a = NSGA2::NSGA2_random_generator.get_RandomInteger(0, this->Population_Size - 1);
		int b = NSGA2::NSGA2_random_generator.get_RandomInteger(0, this->Population_Size - 1);
		Individual* child_a = new Individual(this->Number_Variables, this->Number_ObjectFunction, 0);
		Individual* child_b = new Individual(this->Number_Variables, this->Number_ObjectFunction, 0);
		this->simulated_binary_x(*this->population[a], *this->population[b], *child_a, *child_b);
		this->polynomial_mutation(*child_a);
		this->polynomial_mutation(*child_b);
		this->prerant_chirdern.push_back(child_a);
		this->prerant_chirdern.push_back(child_b);
	}
	for (int i = 0; i < this->Population_Size; i++)
	{
		Individual* ones = new Individual(*this->population[i]);
		this->prerant_chirdern.push_back(ones);
	}
	this->nondominated_sort(this->prerant_chirdern);
	this->calculate_croweddist(this->prerant_chirdern);
	this->population.clear();
	for (int i = 0; i < this->Population_Size; i++)
	{
		Individual* ones = new Individual(*this->prerant_chirdern[i]);
		this->population.push_back(ones);
	}
}

void NSGA2::run(int iteration_times)
{
	for (int i = 0; i < iteration_times; i++)
	{
		printf_s("for iteration times : %d\n", i);
		this->evoluate(0);
	}
}

void NSGA2::display()
{
	ofstream out;
	out.open("log.txt", ios::out);
	for (int i = 0; i < this->Population_Size; i++)
		out << this->population[i]->ObjectFunction_Value[0] << " ";
	out << "\n";
	for (int i = 0; i < this->Population_Size; i++)
		out << this->population[i]->ObjectFunction_Value[1] << " ";
	out << '\n';
	for (int i = 0; i < this->Population_Size; i++)
		out << this->population[i]->ObjectFunction_Value[2] << " ";
	out << "\n";
	out.close();
}

#endif // !__NEGA2__

