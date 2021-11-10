#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <set>
#include <cmath>
#define RIGOR 0.001f

class polinom
{
public:
	polinom();
	polinom(const std::string &str);
	polinom(const polinom &other);
    polinom &operator = (const polinom &other);
	polinom operator + (polinom p);
	polinom operator - (polinom p);
	polinom operator * (polinom p);
	polinom operator / (polinom p);
	polinom operator % (polinom p);
	float from(double x) const;
	polinom derivative (void);
	std::vector<double> 	coeff;
	int 					degree;
};

class Shturm
{
public:
	Shturm(polinom p);
	std::vector<polinom> f;
	void print();
	int from(float x);
	double minus_inf;
	double plus_inf;
	std::vector<double>::iterator if_notcorrect_set();
	void scope();
	std::vector<double> l;
	double find(int i);
};

polinom::polinom()
{
	degree = 0;
}

std::map<int, double> parser_coeff(const std::string &str)
{
	size_t i = 0;
	size_t j = 0;
	size_t g = 0;
	std::string copy;
	std::vector<std::string> vs;
	std::map<int, double> coeff;
	copy = str;
	while(copy.length() > 0)
	{
		i = copy.find("+", 1);
		j = copy.find("-", 1);
		if (i == std::string::npos && j == std::string::npos)
		{
			vs.push_back(copy);
			break ;
		}
		if (i < j)
			g = i;
		else
			g = j;
		vs.push_back(std::string(copy.begin(), copy.begin() + g));
		copy = std::string(copy.begin() + g, copy.end());
	}
	int a;
	double f;
	for (i = 0; i < vs.size(); ++i)
	{
		j = vs[i].find("x^");
		if (j != std::string::npos)
		{
			a = atoi(std::string(vs[i].begin() + j + 2, vs[i].end()).c_str());
			f = atof(std::string(vs[i].begin(), vs[i].begin() + j).c_str());
			if (f == 0.0f && std::string(vs[i].begin(), vs[i].begin() + j).find("0") == std::string::npos)
			{
				if (std::string(vs[i].begin(), vs[i].begin() + j).find("-") != std::string::npos)
					f = -1;
				else
					f = 1;
			}
		}
		else
		{
			j = vs[i].find("x");
			if (j != std::string::npos)
			{
				a = 1;
				f = atof(std::string(vs[i].begin(), vs[i].begin() + j).c_str());
				if (f == 0.0f && std::string(vs[i].begin(), vs[i].begin() + j).find("0") == std::string::npos)
				{
					if (std::string(vs[i].begin(), vs[i].begin() + j).find("-") != std::string::npos)
						f = -1;
					else
						f = 1;
				}
			}
			else
			{
				a = 0;
				f = atof(vs[i].c_str());
			}
		}
		coeff[a] += f;
	}
	return (coeff);
}

std::string parser(const std::string &str)
{
	size_t i = 0;
	std::string copy;
	std::string left;
	std::string right;
	std::string answer;
	std::vector<std::string> vs;
	std::map<int, double> left_coeff;
	std::map<int, double> right_coeff;
	for(i = 0; i < str.length(); ++i)
		if (str[i] != ' ')
			copy += str[i];
	if ((i = copy.find("=")) == std::string::npos)
	{
		std::cout << "error: no equal\n";
		exit(-1);
	}
	else
	{
		left = std::string(copy.begin(), copy.begin() + i);
		right = std::string(copy.begin() + i + 1, copy.end());
		if (left.size() == 0 || right.size() == 0)
		{
			std::cout << "error\n";
			exit(-1);
		}
		left_coeff = parser_coeff(left);
		right_coeff = parser_coeff(right);
	}
	for (std::map<int, double>::iterator it = right_coeff.begin(); it != right_coeff.end(); ++it)
		left_coeff[it->first] -= it->second;

	for (std::map<int, double>::iterator it = left_coeff.begin(); it != left_coeff.end(); ++it)
	{
		if (it->second > 0)
			answer += "+";
		answer += std::to_string(it->second) + "*x^" + std::to_string(it->first);
	}
	return (answer);
}

polinom::polinom(const std::string &str)
{
	size_t i = 0;
	size_t idx = 0;
	std::string k;
	std::string copy;
	copy = str;
	while(copy.length() > 0)
	{
		k = "x^" + std::to_string(idx++);
		if ((i = copy.find(k)) != std::string::npos)
		{
			coeff.push_back(atof(std::string(copy.begin(), copy.begin() + i).c_str()));
			copy = std::string(copy.begin() + i + k.length(), copy.end());
			degree = idx-1;
		}
		else
		{
			coeff.push_back(0.0f);
		}
	}
	while (degree > 0 && coeff[degree] == 0.0)
	{
		coeff.pop_back();
		degree--;
	}
	
}

polinom::polinom(const polinom &other)
{
	*this = other;
}

polinom &polinom::operator = (const polinom &other)
{
	if (this == &other)
		return (*this);
	coeff = other.coeff;
	degree = other.degree;
	return (*this);
}

std::ostream& operator<<(std::ostream &out, const polinom &p)
{
	out << "polinom: ";
	for(int i = 0; i <= p.degree; i++)
	{
		if(i > 0 && p.coeff[i] >= 0)
			out << "+";
		out << p.coeff[i] << "*x^"<< i;
	}
	out << "=0\ndegree: " << p.degree << "\n";
	return (out);
}

float polinom::from(double x) const
{
	double s = 0.0;

	for(int i = degree; i >= 0; i--)
		s = s * x + coeff[i];
	return (s);
}

polinom polinom::operator + (polinom p)
{
    polinom sum;
	
	for (sum.degree = 0; sum.degree <= degree && sum.degree <= p.degree; sum.degree++)
		sum.coeff.push_back(coeff[sum.degree] + p.coeff[sum.degree]);
	for (; sum.degree <= degree || sum.degree <= p.degree; sum.degree++)
	{
		if (sum.degree <= degree)
			sum.coeff.push_back(coeff[sum.degree]);
		else
			sum.coeff.push_back(p.coeff[sum.degree]);
	}
	sum.degree--;
	while (sum.degree > 0 && abs(sum.coeff[sum.degree]) <= RIGOR)
	{
		sum.coeff.pop_back();
		sum.degree--;
	}
	return(sum);
}

polinom polinom::operator - (polinom p)
{
    polinom dif;
	
	for (dif.degree = 0; dif.degree <= degree && dif.degree <= p.degree; dif.degree++)
		dif.coeff.push_back(coeff[dif.degree] - p.coeff[dif.degree]);
	for (; dif.degree <= degree || dif.degree <= p.degree; dif.degree++)
	{
		if (dif.degree <= degree)
			dif.coeff.push_back(coeff[dif.degree]);
		else
			dif.coeff.push_back(-p.coeff[dif.degree]);
	}
	dif.degree--;
	while (dif.degree > 0 && abs(dif.coeff[dif.degree]) <= 10*RIGOR)
	{
		dif.coeff.pop_back();
		dif.degree--;
	}
	return(dif);
}

polinom polinom::operator * (polinom p)
{
    polinom res;

	for (res.degree = 0; res.degree <= degree + p.degree; res.degree++)
		res.coeff.push_back(0);
	res.degree--;
	for (int i = 0; i <= degree; i++)
		for (int j = 0; j <= p.degree; j++)
			res.coeff[i + j] += coeff[i] * p.coeff[j]; 
	return (res);
}

polinom polinom::operator / (polinom p)
{
    polinom res;
	polinom temp = *this;

	if (p.degree > degree)
		return (res);
	for (res.degree = 0; res.degree <= degree - p.degree; res.degree++)
		res.coeff.push_back(0);
	res.degree--;
	int i;
	while ((i = temp.degree) >= p.degree)
	{
		res.coeff[i - p.degree] = temp.coeff[temp.degree] / p.coeff[p.degree];
		temp  = temp - p * polinom(std::to_string(res.coeff[i - p.degree])+"*x^" + std::to_string(i - p.degree));
		if (temp.degree == 0 && abs(temp.coeff[0]) <= RIGOR)
			break ;
	}
	return (res);
}

polinom polinom::operator % (polinom p)
{
	return (*this - (*this / p) * p);
}

polinom polinom::derivative ()
{
	polinom res;

	for (res.degree = 0; res.degree < degree; res.degree++)
		res.coeff.push_back(coeff[res.degree + 1] * (res.degree + 1.0));
	if (res.degree > 0)
		res.degree--;
	return (res);
}

bool  operator == (polinom p1, polinom p2)
{
	return (p1.coeff == p2.coeff && p1.degree == p2.degree);
}

bool  operator != (polinom p1, polinom p2)
{
	return (!(p1 == p2));
}

bool if_null(polinom p)
{
	for (int i = 0; i < p.coeff.size(); ++i)
		if (abs(p.coeff[i]) > RIGOR)
			return (0);
	return (1);
}

polinom nod (polinom p1, polinom p2)
{
	polinom nul;
	polinom del;
	polinom a;

	nul.coeff.push_back(0.0f);
	if (p1 == p2)
		return (p1);
	if (p2.degree > p1.degree)
	{
		a = p1;
		p1 = p2;
		p2 = a;
	}
	while (if_null(p1 % p2) == 0)
	{
		
		for (int i = 0; i < p1.coeff.size(); ++i)
			p1.coeff[i] /= p1.coeff[p1.coeff.size()-1];
		for (int i = 0; i < p2.coeff.size(); ++i)
			p2.coeff[i] /= p2.coeff[p2.coeff.size()-1];
		a = p2;
		p2 = p1 % p2;
		p1 = a;
	}
	return (p2);
}

int Shturm::from(float x)
{
	int	count = 0;
	float a = 0.0f;

	for (int i = 0; i < f.size(); ++i)
	{
		if (abs(f[i].from(x)) < RIGOR)
			continue ;
		if ((f[i].from(x) > RIGOR && a < 0) || (f[i].from(x) < -RIGOR && a > 0))
			count++;
		a = f[i].from(x);
	}
	return(count);
}

Shturm::Shturm(polinom p)
{
	polinom nul;
	polinom temp;

	plus_inf = 1;
	nul.coeff.push_back(0.0f);
	f.push_back(p / nod(p, p.derivative()));
	f.push_back(f[0].derivative());
	while (if_null(temp = (nul - f[f.size()-2] % f[f.size()-1])) == 0)
		f.push_back(temp);
	for (int i = 0; i < f[0].degree; ++i)
		if ((f[0].coeff[f[0].degree-1-i] / f[0].coeff[f[0].degree]) > 1)
			plus_inf += (f[0].coeff[f[0].degree-1-i] / f[0].coeff[f[0].degree]) * (f[0].coeff[f[0].degree-1-i] / f[0].coeff[f[0].degree]);
		else
			plus_inf += f[0].coeff[f[0].degree-1-i] * f[0].coeff[f[0].degree-1-i];
	minus_inf = - plus_inf;

	l.push_back(minus_inf);
	l.push_back(plus_inf);
}

void Shturm::print()
{
	for(int i = 0; i < f.size(); ++i)
		std::cout << f[i];
}


std::vector<double>::iterator Shturm::if_notcorrect_set()
{
	for (std::vector<double>::iterator it = l.begin(); it != l.end() - 1; ++it)
	{
		if (from(*it) - from(*(it + 1)) > 1)
			return (it + 1);
	}
	for (std::vector<double>::iterator it = l.begin(); it != l.end() - 1; )
	{
        if (from(*it) - from(*(it + 1)) == 0)
            l.erase(it + 1);
        else 
        	++it;
    }
	return (l.end());
}

double Shturm::find(int i)
{
	double min_in = l[i];
	double max_in = l[i+1];
	polinom p = f[0];

	double mid_val =  p.from((min_in + max_in) / 2);
	int k = 0;
	int g = 0;
	if (p.from(max_in) > p.from(min_in))
		k = 1;
	while (abs(mid_val) > RIGOR/100000)
	{
		if (mid_val < -RIGOR/100000)
		{
			if (k == 0)
				max_in = (min_in + max_in) / 2;
			else
				min_in = (min_in + max_in) / 2;
		}
		else
		{
			if (k == 0)
				min_in = (min_in + max_in) / 2;
			else
				max_in = (min_in + max_in) / 2;
		}
		mid_val =  p.from((min_in + max_in) / 2);
	}
	return ((min_in + max_in) / 2);
}

void Shturm::scope()
{
	std::vector<double>::iterator it;
	while ((it = if_notcorrect_set()) != l.end())
	{
		l.insert(it, (*it + *(it - 1))/2);
	}
	/*
	for (int i = 0; i < l.size(); ++i)
	{
		std::cout << l[i] << " ";
	}
	std::cout << "\n";
	for (int i = 0; i < l.size(); ++i)
	{
		std::cout << from(l[i]) << " ";
	}
	std::cout << "\n";
	*/
	std::cout << "solution:\n";
	if (l.size() > 1)
		for (int i = 0; i < l.size() - 1; ++i)
			std::cout << "X" << i + 1 << " = " << round(find(i) * 100) / 100 << "\n";
	else
		std::cout << "no real roots of the polynomial\n";
}

double mySqrt(double num) {  
    double x = 1;

    while(abs(x*x - num) >= 0.001 )
        x = ((num/x) + x) / 2;

    return x;

}  



void solution_pol2(polinom p)
{
	double disk = p.coeff[1] * p.coeff[1] - 4 * p.coeff[0] * p.coeff[2];
	std::cout << "Discriminant = " << disk << "\nsolution:\n";
	if (disk < 0)
	{

		double x = -p.coeff[1] / (2.0 * p.coeff[2]);
		double y = mySqrt(-disk) / (2.0 * p.coeff[2]);
		std::cout << "X1 = " << x << " + " << y << "i\n";
		std::cout << "X2 = " << x << " - " << y << "i\n"; 
	}
	else if (disk == 0)
	{
		double x = -p.coeff[1] / (2.0 * p.coeff[0]);
		std::cout << "X1 = X2 = " << x << "\n";
	}
	else
	{
		double x1 = (-p.coeff[1] + mySqrt(disk)) / (2.0 * p.coeff[2]);
		double x2 = (-p.coeff[1] - mySqrt(disk)) / (2.0 * p.coeff[2]);
		std::cout << "X1 = " << x1 << "\n";
		std::cout << "X2 = " << x2 << "\n"; 
	}
}


void solution(polinom p)
{
	std::cout << p;
	if (p.degree == 2)
	{
		solution_pol2(p);
	}
	else if (p.degree > 0)
	{
		Shturm sh(p);
		sh.scope();
	}
	else
	{
		std::cout << "solution:\n";
		if (p.coeff[0] != 0.0f)
			std::cout << "no real roots of the polynomial\n";
		else
			std::cout << "all real numbers are solution\n";
	}
}

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		polinom p1(parser("1 * x^0 - x^1 + 7 * x^2 = 1 * x^1 - 5x + x^3"));
		polinom p2(parser("-2 * x^0 + 1 * x^1 = 0"));
		polinom p3(parser("3*x^0+1*x^2 = 0"));
		polinom p4(parser("-4*x^0+1*x^1 = 0"));
		polinom p5(parser("2*x^0 = 0"));
		polinom p6(parser("42 = 42"));
		solution(p1);
		solution(p2);
		solution(p3);
		solution(p4);
		solution(p5);
		solution(p6);
		solution(p1 * p2 * p3 * p4 * p5);
		solution(p1 * p2 * p3 * p4 * p5 *p1);
	}
	else
	{
		polinom p(parser(std::string(argv[1])));
		solution(p);
	}
	return(0);
}
