#include "pch.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <set>
#include <vector>
#include <map>
#include <climits>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <regex>
#include <random>
#include <queue>
#include <stack>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <cstring>
#include <complex>
#include <cassert>
#include <iterator>
#include <functional>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <chrono>
#include <future>
#include <memory>
#include <csignal>
#include <thread>
#include <condition_variable>
#include "UserInterruptHandler.h"
using namespace std;
using namespace chrono;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);



	return 0;
}
*/
// SCHAB

/*
struct list { int number; list* next; };

int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	list *begin, *temp;
	int x;

	begin = new list; temp = begin;
	cin >> temp->number;

	do
	{
		temp->next = new list;
		temp = temp->next;
		cin >> temp->number;
		temp->next = 0;
	} while (!cin.eof());

	for (temp = begin; temp; temp = temp->next)
		cout << temp->number << " ";

	return 0;
}
*/
// STRUCT 1

/*
//struct list { int number; list* next; };

int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	list *begin, *temp, *tmp;
	int a, b, i;
	begin = new list;
	cin >> b;
	temp = begin;
	while (1)
	{
		a = b;
		cin >> b;
		if (cin.eof())break;

		temp->next = new list;
		temp = temp->next;
		temp->number = a;
		temp->next = 0;
	}
	tmp = begin; begin = begin->next;
	delete tmp;

	for (temp = begin, i = 0; i < b - 1 && temp; i++, temp = temp->next);
	tmp = new list;

	tmp->next = temp->next;
	tmp->number = temp->number;
	temp->number = a;
	temp->next = tmp;

	for (temp = begin; temp; temp = temp->next)
		cout << temp->number << " ";

	return 0;
}
*/
// STRUCT 2

/*
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	int n, m, i, j, **a;

	cin >> n >> m;
	a = new int *[n];

	for (i = 0; i < n; i++)
	{
		a[i] = new int[m];
		for (j = 0; j < m; j++)
			cin >> a[i][j];
	}
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		{
			if (i > 0 && j > 0)
			{
				a[i][j] += a[i - 1][j] < a[i][j - 1] ? a[i - 1][j] : a[i][j - 1];
				continue;
			}
			if (i > 0)
			{
				a[i][j] += a[i - 1][j]; continue;
			}
			if (j > 0)
			{
				a[i][j] += a[i][j - 1];
			}
		}
	cout << a[n - 1][m - 1];

	return 0;
}
*/
// TURTLE

/*
#define Pi 3.1415926
using namespace std;

double length(int x1, int y1, int x2, int y2)
{
	return sqrt((double)(x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
}


int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	double l[3], alpha, r, cos, t;
	int i, x[3], y[3];

	for (i = 0; i < 3; i++)
		cin >> x[i] >> y[i];

	l[0] = length(x[1], y[1], x[0], y[0]);
	l[1] = length(x[2], y[2], x[1], y[1]);
	l[2] = length(x[2], y[2], x[0], y[0]);

	cos = (l[0] * l[0] + l[1] * l[1] - l[2] * l[2]) / (2 * l[0] * l[1]);
	alpha = acos(cos);
	t = 180 / ((Pi) / alpha);

	cos = (l[1] * l[1] + l[2] * l[2] - l[0] * l[0]) / (2 * l[1] * l[2]);
	alpha = acos(cos);
	r = 180 / ((Pi) / alpha);

	if (r > t)
		t = r;

	cos = (l[0] * l[0] + l[2] * l[2] - l[1] * l[1]) / (2 * l[0] * l[2]);
	alpha = acos(cos);
	r = 180 / ((Pi) / alpha);

	if (r > t)
		t = r;

	cout.precision(6);

	cout << t;

	return 0;
}
*/
// GEOMETRY 1

/*
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	int k, W, *w, **a, n, s, i, j;
	cin >> k >> W;
	w = new int[k + 1];

	w[0] = 0;
	for (i = 1; i <= k; i++)
		cin >> w[i];

	a = new int*[k + 1];

	for (i = 0; i <= k; i++)
		a[i] = new int[W + 1];

	for (n = 0; n <= W; n++)
		a[0][n] = 0;

	for (n = 1; n <= k; n++)
		a[n][0] = 0;

	for (s = 1; s <= k; s++)
	{
		for (n = 0; n <= W; n++)
		{
			a[s][n] = a[s - 1][n];
			if (n >= w[s] && (a[s - 1][n - w[s]] + w[s] > a[s][n]))a[s][n] = a[s - 1][n - w[s]] + w[s];
		}

		if (a[s][W] == W)
		{
			k = s; break;
		}
	}
	cout << a[k][W];
}
*/
// BACKPACK 1

/*
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	int k, W, *w, *p, **a, n, s, i, j;

	cin >> k >> W;
	p = new int[k + 1];

	p[0] = 0;

	for (i = 1; i <= k; i++)
		cin >> p[i];

	w = new int[k + 1];

	w[0] = 0;

	for (i = 1; i <= k; i++)
		cin >> w[i];

	a = new int*[k + 1];
	for (i = 0; i <= k; i++)
		a[i] = new int[W + 1];

	for (n = 0; n <= W; ++n)
		a[0][n] = 0;


	for (n = 1; n <= k; ++n)
		a[n][0] = 0;

	for (s = 1; s <= k; ++s)
	{
		for (n = 0; n <= W; ++n)
		{
			a[s][n] = a[s - 1][n];
			if (n >= w[s] && (a[s - 1][n - w[s]] + p[s] > a[s][n]))
			{
				a[s][n] = a[s - 1][n - w[s]] + p[s];

			}
		}

	}
	cout << a[k][W];

	return 0;
}
*/
// BACKPACK 2

/*
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);


	int m, n, i, j, **a, **b;
	cin >> n >> m;
	a = new int *[n];
	b = new int *[n + 1];


	for (i = 0; i < n; i++)
	{
		a[i] = new int[m];

		for (j = 0; j < m; j++)
			cin >> a[i][j];
	}

	for (i = 0; i <= n; i++)
	{
		b[i] = new int[m + 1];
		for (int j = 0; j <= m; j++)
			b[i][j] = 0;
	}

	b[1][1] = 1;

	for (i = 2; i <= n; i++)
		for (j = 2; j <= m; j++)
			b[i][j] = b[i - 1][j - 2] + b[i - 2][j - 1];

	if (!b[n][m])
	{
		cout << "-";
		return 0;
	}

	for (i = 1; i < n; i++)
		for (j = 1; j < m; j++)


	{
		if (i - 1 == 0 && j - 2 > 0)
			continue;
		if (i - 1 >= 0 && j - 2 >= 0 && i - 2 >= 0 && j - 1 >= 0)
		{
			a[i][j] = a[i][j] + (a[i - 1][j - 2] >= a[i - 2][j - 1] ? a[i - 1][j - 2] : a[i - 2][j - 1]);
			continue;
		}

		if (i - 1 >= 0 && j - 2 >= 0)
		{
			a[i][j] = a[i][j] + a[i - 1][j - 2];
			continue;
		}

		if (i - 2 >= 0 && j - 1 >= 0)
		{
			a[i][j] = a[i][j] + a[i - 2][j - 1];
		}
	}

	cout << a[n - 1][m - 1];

	return 0;
}
*/
// HORSE

/*
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	int  n, i, j;
	cin >> n;
	int *a = new int[n + 1];

	a[0] = 0;

	for (int i = 1; i <= n; i++)
	{
		cin >> a[i];
	}

	for (i = 3; i <= n; i++)
	{
		a[i] = a[i] + (a[i - 2] > a[i - 3] ? a[i - 2] : a[i - 3]);
	}
	cout << a[n];

	return 0;
}
*/
// GRASSHOPPER

/*
#define _CRT_SECURE_NO_DEPRECATE
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <set>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>
#include <random>
#include <queue>
#include <stack>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <cstring>
#include <complex>
#include <cassert>
#include <iterator>
#include <functional>
typedef long long int64;
typedef long double ld;
const int inf = (int)2e+9 + 2;
const int mod = (int)1e+9 + 7;
const int h_prime = 1103;
const double eps = 1e-6;
//const int maxn = 1e+6;


#define pb push_back
using namespace std;
typedef pair<ld, ld> vec;
typedef vector <int64> vint;
typedef vector <vint> vvint;
vint primes;
vint HH, pp;

ld dist(vec a, vec b)
{
	return sqrtl((a.first - b.first)*(a.first - b.first) + (a.second - b.second)*(a.second - b.second));
}
ld len(vec x)
{
	return sqrtl(x.first*x.first + x.second*x.second);
}
ld angle(vec a, vec b)
{
	return acos((a.first * b.first + a.second*b.second) / (len(a)*len(b)));
}
ld ar_tr(vec a, vec b)
{
	return a.first*b.second - b.first*a.second;
}

int64 get_hash(int l, int r)
{
	int64 res = HH[r];
	if (l)
		res -= HH[l - 1] * pp[r - l + 1];
	return res;
}

void init_hashes(string s)
{
	int64 hash = s[0];
	HH.pb(hash);
	int64 ppow = 1;
	for (size_t i = 0; i < s.length(); ++i)
	{
		pp.pb(ppow);
		ppow *= h_prime;
	}
	for (size_t i = 1; i < s.length(); ++i)
	{
		hash *= h_prime;
		hash += s[i];
		HH.pb(hash);
	}
}

inline int64 gcd(int64 x, int64 y)
{
	while (y)
	{
		x %= y;
		swap(x, y);
	}
	return x;
}
void init_primes(int64 maxn)
{
	vector <bool> np(maxn, false);
	np[0] = np[1] = true;
	for (int i = 2; i*i < maxn; ++i)
		if (!np[i])
		{
			for (int j = i * i; j < maxn; j += i)
				np[j] = true;
			primes.pb(i);
		}
}
int64 binpow(int64 a, int64 b)
{
	int64 res = 1;
	while (b)
	{
		if (b & 1)
			res *= a;
		a *= a;
		b /= 2;
	}
	return res;
}

#ifdef TREE
vector<int64> tree, nval;

void push(int64 v) {
	if (nval[v]) {
		tree[v << 1] += nval[v];
		tree[v << 1 | 1] += nval[v];
		nval[v << 1] += nval[v];
		nval[v << 1 | 1] += nval[v];
		nval[v] = 0;
	}
}

void upd(int64 v, int64 tl, int64 tr, int64 l, int64 r, int64 val) {
	int64 tm = (tl + tr) >> 1;
	push(v);
	if (l > r) return;
	if (l == tl && r == tr) {
		tree[v] += val;
		nval[v] += val;
		return;
	}
	upd(v << 1, tl, tm, l, min(r, tm), val);
	upd(v << 1 | 1, tm + 1, tr, max(l, tm + 1), r, val);
	tree[v] = tree[v << 1] + tree[v << 1 | 1];
}

int64 get(int64 v, int64 tl, int64 tr, int64 l, int64 r) {
	push(v);
	int64 tm = (tl + tr) >> 1;
	if (l > r)
		return int64(1e9);
	if (l == tl && r == tr)
		return tree[v];
	return get(v << 1, tl, tm, l, min(r, tm) + get(v << 1 | 1, tm + 1, tr, max(l, tm + 1), r);
}
#endif

struct pt {
	double x, y;
	pt(int a, int b)
	{
		x = a;
		y = b;
	}
};

bool cmp(pt a, pt b) {
	return a.x < b.x || a.x == b.x && a.y < b.y;
}

bool cw(pt a, pt b, pt c) {
	return a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y) < 0;
}

bool ccw(pt a, pt b, pt c) {
	return a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y) > 0;
}

void convex_hull(vector<pt> & a) {
	if (a.size() == 1)  return;
	sort(a.begin(), a.end(), &cmp);
	pt p1 = a[0], p2 = a.back();
	vector<pt> up, down;
	up.push_back(p1);
	down.push_back(p1);
	for (size_t i = 1; i < a.size(); ++i) {
		if (i == a.size() - 1 || cw(p1, a[i], p2)) {
			while (up.size() >= 2 && !cw(up[up.size() - 2], up[up.size() - 1], a[i]))
				up.pop_back();
			up.push_back(a[i]);
		}
		if (i == a.size() - 1 || ccw(p1, a[i], p2)) {
			while (down.size() >= 2 && !ccw(down[down.size() - 2], down[down.size() - 1], a[i]))
				down.pop_back();
			down.push_back(a[i]);
		}
	}
	a.clear();
	for (size_t i = 0; i < up.size(); ++i)
		a.push_back(up[i]);
	for (size_t i = down.size() - 2; i > 0; --i)
		a.push_back(down[i]);
}
bool line(double A, double B, double C, double x, double y)
{
	return (A*x + B * y + C > 0);
}
bool check(double A, double B, double C, vector <pt>& a, vector <pt>& b)
{
	bool f = line(A, B, C, a[0].x, a[0].y);
	for (int i = 0; i < a.size(); ++i)
	{
		if (line(A, B, C, a[i].x, a[i].y) != f)
			return false;
	}
	for (int i = 0; i < b.size(); ++i)
	{
		if (line(A, B, C, b[i].x, b[i].y) == f)
			return false;
	}
	return true;
}
int main()
{
	ios::sync_with_stdio(0);
	cin.tie(0), cout.tie(0);
	//#ifdef _DEBUG
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	//#endif
	vector <pt> p, b;
	int n;
	cin >> n;
	pair <int, int> mina = { inf,inf }, minb{ inf,inf }, maxa = { -inf,-inf }, maxb = { -inf,-inf };
	for (int i = 0; i < n; ++i)
	{
		int x, y, t;
		cin >> x >> y >> t;
		if (t == 1)
		{
			mina.first = min(mina.first, x);
			mina.second = min(mina.second, y);
			maxa.first = max(maxa.first, x);
			maxa.second = max(maxa.second, y);
			p.pb({ x,y });
		}
		else
		{

			minb.first = min(minb.first, x);
			minb.second = min(minb.second, y);
			maxb.first = max(maxb.first, x);
			maxb.second = max(maxb.second, y);
			b.pb({ x,y });
		}
	}
	if (maxa.first < minb.first || maxb.first < mina.first || maxa.second < minb.second || maxb.second < mina.second)
	{
		cout << "Yes";
		return 0;
	}
	convex_hull(p);
	convex_hull(b);
	for (int i = 0; i < p.size(); ++i)
	{
		pt v = p[i];
		pt next = p[(i + 1) % p.size()];
		ld A = v.y - next.y;
		ld B = next.x - v.x;
		ld C = v.x*next.y - next.x*v.y;
		if (check(A, B, C, p, b))
		{
			cout << "Yes";
			return 0;
		}
	}
	for (int i = 0; i < b.size(); ++i)
	{
		pt v = b[i];
		pt next = b[(i + 1) % b.size()];
		ld A = v.y - next.y;
		ld B = next.x - v.x;
		ld C = v.x*next.y - next.x*v.y;
		if (check(A, B, C, b, p))
		{
			cout << "Yes";
			return 0;
		}
	}
	if (p.size() == 1 && b.size() == 1)
	{
		cout << "Yes";
		return 0;
	}
	cout << "No";
}
*/
// GG

/*
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);


	int **a, n, m, i, j, x, y;

	cin >> n >> m;

	a = new int*[n];

	for (i = 0; i < n; i++)
	{
		a[i] = new int[n];
		for (j = 0; j < n; j++)
		{
			a[i][j] = 0;
		}
	}

	for (i = 0; i < m; i++)
	{
		cin >> x >> y;

		a[x - 1][y - 1] = 1;
		a[y - 1][x - 1] = 1;
	}

	for (i = 0; i < n; i++)
	{
		cout << i + 1 << ": ";
		for (j = 0; j < n; j++)
		{
			if (a[i][j] == 1)
			{
				cout << j + 1 << " ";

			}
		}
		cout << endl;
	}

	return 0;
}
*/
// GRAPH

/*
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	int **a,*v0,*v1, n, m, i, j, x, y;

	cin >> n >> m;

	a = new int*[n];

	for (i = 0; i < n; i++)
	{
		a[i] = new int[n];
		for (j = 0; j < n; j++)
		{
			cin >> a[i][j];
		}
	}

	v0 = new int[n];
	v1 = new int[n];

	for (i = 0; i < n; i++)
	{
		v1[i] = a[m - 1][i];
	}

	while (1)
	{
		for (i = 0; i < n; i++)
		{
			v0[i] = v1[i];
		}

		for (i = 0; i < n; i++)
		{
			if (v0[i] == 1)
			{
				for (j = 0; j < n; j++)
				{
					if (v1[j] == 0)
					{
						v1[j] = a[i][j];
					}
				}
			}
		}

		for (i = 0; i < n && (v1[i] == v0[i]); i++);

		if (i == n)
		{
			break;
		}
	}

	cout << m << " ";

	for (i = 0; i < n; i++)
	{
		if ((i+1)!=m && v1[i] == 1)
		{
			cout << i + 1 << " ";
		}
	}

	return 0;
}
*/
// GRAPH 2

/*
const int INF = 1000000000;
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	int n, s, f;
	cin >> n >> s >> f;

	vector < vector < pair<int, int> > > g(n);



	vector<int> d(n, INF), p(n);
	d[s] = 0;
	vector<char> u(n);
	for (int i = 0; i < n; ++i)
	{

		int v = -1;
		for (int j = 0; j < n; ++j)
			if (!u[j] && (v == -1 || d[j] < d[v]))
				v = j;
		if (d[v] == INF)
			break;
		u[v] = true;

		for (size_t j = 0; j < g[v].size(); ++j)
		{

			int to = g[v][j].first,
				len = g[v][j].second;
			if (d[v] + len < d[to])
				{
				d[to] = d[v] + len;
				p[to] = v;
			}
		}
	}
	copy(d.begin(), d.end(), ostream_iterator<int>(cout, " "));

	return 0;
}
*/
// ???

/*
int NextPerest(int *a, int n)
{
	int i, j, k;
	for (i = n - 2; i >= 0; i--)
		if (a[i] < a[i + 1])
			break;

	if (i < 0)
		return -1;
	for (j = n - 1; a[j] < a[i]; j--);

	k = a[i];
	a[i] = a[j];
	a[j] = k;
	for (j = i + 1, i = n - 1; j < i; j++, i--)
	{
		k = a[j];
		a[j] = a[i];
		a[i] = k;
	}

	return 1;
}

int cost(int **c, int *p, int n)
{
	int s = 0, i;
	for (i = 0; i < n - 1; i++)
		s = s + c[p[i] - 1][p[i + 1] - 1];
	return s;
}

int main()
{

	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	int i, **c, n, j;
	int *p, *pm, sm, s, r, z;
	cin >> n;
	c = new int*[n];

	for (i = 0; i < n; i++)
	{
		c[i] = new int[n];
		for (j = 0; j < n; j++)
			cin >> c[i][j];
	}

	p = new int[n];
	pm = new int[n];
	for (i = 0; i < n; i++)
		p[i] = i + 1;

	sm = cost(c, p, n);

	for (j = 0; j < n; j++)
	{
		pm[j] = p[j];
	}

	while (1)
	{
		if (NextPerest(p, n) == -1)
		{
			break;
		}

		s = cost(c, p, n);
		if (s < sm)
		{
			sm = s;
			for (j = 0; j < n; j++)
				pm[j] = p[j];
		}
	}
	cout << sm << endl;
	for (j = 0; j < n; j++)
		cout << pm[j] << " ";

	return 0;

}
*/
// KAMEN ZHORA

/*
int *c;
int *path;
int **a;

int gamilton(int **a, int n, int k, int v0)
{
	int v, q1 = 0;
	for (v = 0; v < n && !q1; v++)
	{
		if (a[v][path[k - 1]] || a[path[k - 1]][v])
		{
			if (k == n && v == v0)
				q1 = 1;
			else if (c[v] == -1)
			{
				c[v] = k;
				path[k] = v;
				q1 = gamilton(a, n, k + 1, v0);
				if (!q1) c[v] = -1;
			}
			else continue;
		}
	}
	return q1;
}

int main()
{

	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	int n, m, i, j, v0;

	cin >> n >> m;

	v0 = m - 1;
	a = new int *[n];

	for (i = 0; i < n; i++)
	{
		a[i] = new int[n];
		for (j = 0; j < n; j++)
			cin >> a[i][j];
	}

	c = new int[n];

	for (j = 0; j < n; j++)
		c[j] = -1;

	path = new int[n];
	path[0] = v0;
	c[v0] = v0;

	if (gamilton(a, n, 1, v0))
	{
		for (i = 0; i < n; i++)
			cout << path[i] + 1 << " ";

		cout << path[0] + 1;
	}
	else
	{
		cout << -1;
	}

	return 0;
}
*/
// N 211

/*
int main()
{

	const int maxInt = 0x7fffffff;


	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	int n, m, i, j, **a, k, l;

	cin >> n >> m >> l;

	a = new int *[n];

	for (i = 0; i < n; i++)
	{
		a[i] = new int[n];
	}


	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			cin >> a[i][j];
			if (a[i][j] == -1)
			{
				a[i][j] = maxInt;
			}
		}


	for (k = 0; k < n; k++)
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				if (a[i][k] < maxInt && a[k][j] < maxInt)
				{
					a[i][j] = (a[i][j] < a[i][k] + a[k][j]) ? a[i][j] : (a[i][k] + a[k][j]);
				}

	if (a[m - 1][l - 1] == maxInt)
	{

		cout << -1;
	}
	else
	{

		cout << a[m - 1][l - 1] << " ";
	}

	return 0;

}
*/
// FLOID 

/*
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	const int Max_Int = 0x7fffffff;
	int n, s, f, v;

	cin >> n >> s >> f;

	int *a, *l;

	a = new int[n];
	l = new int[n];

	for (int i = 0; i < n; i++)
	{
		a[i] = 0;
		l[i] = Max_Int;
	}

	int **w = new int*[n];

	for (int i = 0; i < n; i++)
	{
		w[i] = new int[n];
		for (int j = 0; j < n; j++)
		{
			cin >> w[i][j];
		}
	}

	l[s - 1] = 0;

	for (int i = 0; i < n; i++)
	{
		v = -1;

		for (int j = 0; j < n; j++)
		{
			if (!a[j] && (v == -1 || l[j] < l[v]))
			{
				v = j;
			}
		}

		if (l[v] == Max_Int)
		{
			break;
		}

		a[v] = 1;

		for (int j = 0; j < n; j++)
		{
			if ((w[v][j] != -1) && (l[v] + w[v][j] < l[j]))
			{
				l[j] = l[v] + w[v][j];
			}
		}
	}

	cout << l[f - 1];

	return 0;
}
*/
// DXTRE

/*
struct elemSt { elemSt*prev; char c; };
elemSt * Push(elemSt *b, char c)
{
	elemSt *t = new elemSt;
	t->prev = b;
	t->c = c;
	return t;
}

elemSt* Pop(elemSt *b)
{
	if (!b)
		return 0;
	elemSt *t = b;
	b = b->prev;
	delete t;
	return b;
}

char Top(elemSt *b)
{
	if (!b)
		return 0;
	return b->c;
}

void Out(elemSt *b)
{
	while (b)
	{
		cout << b->c << ' ';
		b = b->prev;
	}
}

elemSt*freeMem(elemSt*begin)
{
	elemSt *t;
	while (begin)
	{
		t = begin; begin = begin->prev;
		delete t;
	}
	return begin;
}
int bracket(char c)
{
	switch (c)
	{
	case '(':return 1;
	case '{':return 2;
	case '[':return 3;
	case '<':return 4;
	case ')':return 5;
	case '}':return 6;
	case ']':return 7;
	case '>':return 8;
	}

	return 0;
}

int main()
{

	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	string s;
	char si, sst;
	int i, k, l;
	elemSt *begin = 0;

	cin >> s;

	//while (si = cin.get(), !cin.eof())
	for (i = 0; i < s.length(); i++)
	{
		si = s[i];
		l = bracket(si);
		if (l <= 4 && l > 0)
		{
			begin = Push(begin, si);
			continue;
		}
		if (!begin)
		{
			cout << 0;
			return 0;
		}

		sst = Top(begin);
		k = bracket(sst);

		if (k + 4 != l)
		{
			cout << 0;
			begin = freeMem(begin);
			return 0;
		}
		begin = Pop(begin);
	}

	if (!begin)
		cout << 1;
	else
	{
		cout << 0;
		begin = freeMem(begin);
	}
	return 0;
}
*/
// STACK

/*
int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	string n, z;

	cmatch result;
	regex a("([a-zA-Z0-9]+[\\s]*=[\\s]*[a-zA-Z0-9]+;)");
	while (getline(cin, n))
	{
		bool k = true;
		while (k)
		{
			if (regex_search(n.c_str(), result, a))
			{
				cout << result.str() << endl;
				int i = n.find(result.str());
				string m = result.str();
				n.erase(i, m.length());
			}
			else
				k = false;
		}


	}
}
*/
// FIRST 

/*
int input()
{
	int n;
	cin >> n;
	return n;
}

int summerize(int var, int sum)
{
	return var + sum;
}

int func(int n, int sum, int index)
{
	if (n != index)
	{
		return func(n, summerize(input(), sum), index + 1);
	}

	return sum;
}

void exec(int sum)
{
	cout << sum;
}


int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	exec(func(input(), 0, 0));
}
*/
// JAJA

/*
double len(int x1, int y1, int x2, int y2)
{
	return sqrt((double)(x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
}

int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	double l[3], xr, yr; int i, x[3], y[3];
	for (i = 0; i < 3; i++)
		cin >> x[i] >> y[i];

	l[0] = len(x[1], y[1], x[0], y[0]);
	l[1] = len(x[2], y[2], x[1], y[1]);
	l[2] = len(x[2], y[2], x[0], y[0]);

	xr = (l[1] * x[0] + l[2] * x[1] + l[0] * x[2]) / (l[1] + l[2] + l[0]);
	yr = (l[1] * y[0] + l[2] * y[1] + l[0] * y[2]) / (l[1] + l[2] + l[0]);

	cout.precision(6);

	cout << xr << " " << yr;

	return 0;
}
*/
// GEO

/*
struct coordinates
{
	int x1, y1, x2, y2, x3, y3;
};

int get_number()
{
	int n;
	cin >> n;
	return n;
}

coordinates get_triangle()
{
	coordinates tr;
	cin >> tr.x1 >> tr.y1 >> tr.x2 >> tr.y2 >> tr.x3 >> tr.y3;
	return tr;
}

vector<coordinates> create_coord_array(int n) {
	vector<coordinates> coords(n);
	return coords;
}

vector<coordinates> triangle_loop(int index, vector<coordinates> triangles)
{
	if (index >= triangles.size())
	{
		return triangles;
	}
	triangles[index] = get_triangle();
	return triangle_loop(index + 1, triangles);

}

coordinates triangle_multiplication(int mul_number, coordinates triangle)
{
	// one particular triangle
	triangle.x1 *= mul_number;
	triangle.y1 *= mul_number;
	triangle.x2 *= mul_number;
	triangle.y2 *= mul_number;
	triangle.x3 *= mul_number;
	triangle.y3 *= mul_number;
	return triangle;
}

vector<coordinates> get_triangles_multiplication(int mul_number, vector<coordinates> triangles, int index) {
	// adding all triangles together
	if (index >= triangles.size()) {
		return triangles;
	}
	triangles[index] = triangle_multiplication(mul_number, triangles[index]);
	return get_triangles_multiplication(mul_number, triangles, index + 1);
}

coordinates count_new_vector(vector<coordinates> triangles, int index) {
	coordinates new_vector;
	new_vector.x1 = triangles[index].x2 - triangles[index].x1;
	new_vector.y1 = triangles[index].y2 - triangles[index].y1;
	// first math vector

	new_vector.x2 = triangles[index].x3 - triangles[index].x2;
	new_vector.y2 = triangles[index].y3 - triangles[index].y2;
	// second math vector

	new_vector.x3 = triangles[index].x1 - triangles[index].x3;
	new_vector.y3 = triangles[index].y1 - triangles[index].y3;
	// third math vector


	return new_vector;
}

vector<coordinates> get_vectors(vector<coordinates> triangles, int index, vector<coordinates> vectors)
{
	if (index >= triangles.size()) {
		return vectors;
	}
	vectors.push_back(count_new_vector(triangles, index)); // add math vectors (of one triangle) to vectors array
	return get_vectors(triangles, index + 1, vectors);

}

void show(vector<coordinates> result, int i)
{
	if (i < result.size())
	{
		cout << result[i].x1 << " " << result[i].y1 << " " << result[i].x2 << " " << result[i].y2 << " " << result[i].x3 << " " << result[i].y3 << endl;
		show(result, i + 1);
	}
}

vector<coordinates> create_empty_vector() {
	vector<coordinates> vectors;
	return vectors;
}

int main()
{
	show(get_vectors(get_triangles_multiplication(get_number(), triangle_loop(0, create_coord_array(get_number())), 0), 0, create_empty_vector()), 0);
}*/
// First 

/*
int get_var()
{
	int n;
	cin >> n;
	return n;
}
vector<int> sortirovka(vector<int> v1, int i, int n, int j, bool indicator)
{
	if (i != n)
	{
		if (n - i - 1 > j)
		{
			if (indicator == true)
			{
				if (v1[j] > v1[j + 1])
					swap(v1[j], v1[j + 1]);
			}
			else {
				if (v1[j] < v1[j + 1])
					swap(v1[j], v1[j + 1]);
			}
			return sortirovka(v1, i, n, j + 1, indicator);
		}
		j = 0;
		return sortirovka(v1, i + 1, n, j, indicator);
	}
	return v1;

}
vector<int> zapolnit_vector(vector<int> v1, int n, int index) {
	if (index != n)
	{
		v1.push_back(get_var());
		return zapolnit_vector(v1, n, index + 1);
	}
	else
		return v1;

}
void vivod(vector<int> v1, vector<int> v2, int n, int i)
{
	if (n > i)
	{
		cout << v1[i] << " " << v2[i] << endl;
		vivod(v1, v2, n, i + 1);
	}
}


int main()
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	int n = get_var();
	vivod(sortirovka(zapolnit_vector(vector<int>(), n, 0), 0, n, 0, true), sortirovka(zapolnit_vector(vector<int>(), n, 0), 0, n, 0, false), n, 0);
}
*/
// Second

/*
int get_var()
{
	int n;
	cin >> n;
	return n;
}

int get_result(int old_number, int new_number, int index)
{
	return new_number * pow(10, index) + old_number;
}

int check_odd_number(int n)
{
	if (n % 2 == 1)
		return n;
	else
		return 0;
}

int get_number(int n)
{
	return n % 10;
}

int func(int n, int result, int index)
{
	if (n == 0)
		return result;
	int res = check_odd_number(get_number(n));
	if (res != 0)
	{
		result = get_result(result, res, index);
		index++;
	}
	return func(n / 10, result, index);
}

void show(int result)
{
	cout << result;
}


int main()
{
	show(func(get_var(),0,0));
}
*/
// Third

/*
int get_var()
{
	int n;
	cin >> n;
	return n;
}
int summarize(int variable, int sum)
{
	return variable + sum;
}
int body(int n, int sum, int index)
{
	if (n != index)
	{
		int a = summarize(get_var(), sum);
		return body(n, a, index + 1);
	}
	else
		return sum;
}

void show(int res)
{
	cout << res;
}

int main()
{
	show(body(get_var(), 0, 0));
}*/
// Fourth

/*mutex mtx;
condition_variable cv; // To avoid a random wake up of Thread.
bool Thatboolean = false;

void threadFunc(shared_ptr<bool> flag, int id) // 
{
	while (!(*flag)) 
	{
		unique_lock<mutex> lk(mtx);  // lock to not come in the Thred func. Not allow to get to the Thred func.
		cv.wait(lk, []() {return Thatboolean; }); // waiting to wake up Thread.
		cout << "Thread Number " + to_string(id) + "\n"; // waiting to get bool true

		if (!*flag)
			Thatboolean = false;
	}
	cout << "Thread " + to_string(id) + " exit" + "\n"; // thread exit.
}

int main() 
{
	UserInterruptHandler h;

	shared_ptr<bool> f = make_shared<bool>(false);
	vector<thread> threads;

	int n;
	cin >> n;

	for (int i = 0; i < n; i++) 
	{
		threads.push_back(thread(threadFunc, f, i));  // for Thread activating function.
	}

	try // waiting 4 error
	{
		while (1) 
		{
			this_thread::sleep_for(chrono::seconds(1)); // break timer for 1 sec.
			Thatboolean = true; // Steins gate
			cv.notify_one();  // wake up Thred.
			h.checkInterruptionAndThrow();
		}
	}

	catch (runtime_error ex) 
	{
		cout << ex.what();
		*f = true;
	}

	Thatboolean = true;
	cv.notify_all(); // all waking up for output.

	for (auto& th : threads) 
	{
		if (th.joinable())
			th.join();
	}

	return 0;
}*/
// FIRST THREAD

/*int n = 16;
mutex mtx;

struct values
{
	int min;    // struct min max value 
	int max;
};

void threadFunc(promise<values> p, vector<long> v1, long  str, long  end)   //
{
	unique_lock<mutex> lock(mtx);   // to separate threads || let thread not come into func.
	values s;

	s.min = INT_MAX;
	s.max = INT_MIN;

	for (long i = str; i < end; i++)   // searching for min and max
	{
		if (v1[i] < s.min) 
		{
			s.min = v1[i];       // finding the min
		}

		if (v1[i] > s.max)
		{
			s.max = v1[i];       // finding the max
		}
	}
	p.set_value(s);     // setting value of promise 
}

int main() 
{
	long int num;

	vector<long> variables;
	cin >> num;

	random_device rd;                                          //
	mt19937 random_engine(rd());                               // RANDOMING THE VALUE FOR THE ARRAY 
	uniform_int_distribution<int> random(INT_MIN, INT_MAX);    //

	for (int i = 0; i < num; i++) {
		variables.push_back(random(random_engine));   // and this too.
	}

	vector<thread> myThread;

	for (int i = 0; i < 5; i++) 
	{
		vector<promise<values>> p(n);    // creating array for promise
		vector<future<values>> f(n);     // creating array for future

		for (int j = 0; j < p.size(); j++) 
		{
			f[j] = p[j].get_future();        // adding promise to future.
		}

		cout << n << ": ";
		auto timeS = high_resolution_clock::now();   // starting timer

		for (int j = 0; j < n; j++)
		{
			myThread.push_back(thread(threadFunc, move(p[j]), variables, j * (num / n), (j + 1) * (num / n)));   // moving(creating) calculating thread line and applying it to the array
		}

		for (int j = 0; j < n; j++)
		{
			values fum = f[j].get();    //  getting value from future 
		}
		auto timeE = high_resolution_clock::now();   // ending timer
		cout << duration_cast<milliseconds>(timeE - timeS).count() << endl;    // output with time used
		n /= 2;


	}

	for (auto& th : myThread)
	{
		if (th.joinable())     // output thread.
			th.join();
	}

	cin >> n;
	return 0;

	//
}*/
// SECOND THREAD

/*condition_variable cv;
condition_variable cvGenerator;

mutex global_mtx;
random_device rd;
mt19937 random_engine(rd());

uniform_int_distribution<int> generator_sleep(1, 3);
uniform_int_distribution<int> random_sleep(10, 15);

bool Thisboolean = true;

struct request
{
	int priority;
	int klass;
};

vector<request> q;
queue<request> temp;

void generateRequest(shared_ptr<bool> flag, int classAm, int qSize) 
{

	mutex mtx;
	unique_lock<mutex> lk(mtx);

	while (!(*flag))
	{

		if (q.size() == qSize) // if queue is Full cout.
		{
			Thisboolean = false;
			cout << "FULL \n";
		}

		cvGenerator.wait(lk, []() { return Thisboolean; }); 
		this_thread::sleep_for(chrono::seconds(generator_sleep(random_engine)));
		uniform_int_distribution<int> random_class(1, classAm);
		uniform_int_distribution<int> random_priority(1, 3);

		request curRequest;

		curRequest.klass = random_class(random_engine);
		curRequest.priority = random_priority(random_engine);

		global_mtx.lock();
		q.push_back(curRequest);
		cv.notify_all();
		global_mtx.unlock();
		cout << "\n \n	? New request created. Group: " + to_string(curRequest.klass) + " priority: " + to_string(curRequest.priority) + "\n";
	}
	cout << "Generator thread exit \n";
}

void deviceFunc(shared_ptr<bool> flag, int id, int group) // 
{
	mutex mtx;
	unique_lock<mutex> lk(mtx);

	bool localCheck = false;
	int r = 0;
	int j;

	while (!(*flag)) 
	{
		if (localCheck == true) 
		{
			cout << "	Thread " + to_string(id + 1) + " Speeping for " + to_string(r) + " seconds \n";

			while (r != 0)
			{
				this_thread::sleep_for(chrono::seconds(1));
				r--;
			}
			cout << "	Thread " + to_string(id + 1) + " Woke up \n";

		}
		localCheck = false;
		j = 0;

		global_mtx.lock();

		for (int i = 0; i < q.size(); i++) // finds specific req with highest priority
		{
			if (q[i].klass == group && q[i].priority > q[j].priority) {
				j = i;
			}
		}

		if (q.size() != 0) // erase suitable request
		{
			q.erase(q.begin() + j);
			localCheck = true;
			Thisboolean = true;
			cvGenerator.notify_one(); // there is free space in queue
			cout << "	! Thread " + to_string(id + 1) + "; Group " + to_string(group) + "; get task: " + to_string(j + 1) + "\n";
		}

		global_mtx.unlock();
		r = random_sleep(random_engine);
	}
	cout << "Thread " + to_string(id) + " exit" + "\n";
}

int main() {
	UserInterruptHandler h;

	shared_ptr<bool> f = make_shared<bool>(false);

	int devicesAmount;
	int groupAmount; 
	int queueSize;

	cin >> queueSize >> groupAmount >> devicesAmount;

	vector<vector<thread>> groups;

	thread generatorThread = thread(generateRequest, f, groupAmount, queueSize);
	for (int i = 0; i < groupAmount; i++) 
	{
		groups.push_back(vector<thread>(devicesAmount));
		for (int j = 0; j < devicesAmount; j++) 
		{
			groups[i][j] = move(thread(deviceFunc, f, (i + 1) * j, i + 1));
		}
	}

	try 
	{
		while (1)
		{
			h.checkInterruptionAndThrow();
		}
	}
	catch (runtime_error ex) 
	{
		cout << ex.what();
		Thisboolean = false;
		*f = true;
	}

	cvGenerator.notify_all();
	cv.notify_all();

	for (int i = 0; i < groupAmount; i++) 
	{
		for (int j = 0; j < devicesAmount; j++) 
		{
			if (groups[i][j].joinable())
				groups[i][j].join();
		}
	}

	if (generatorThread.joinable()) 
	{
		generatorThread.join();
	}
	return 0;
}*/
// THIRD THREAD