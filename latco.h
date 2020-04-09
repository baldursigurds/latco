#ifndef __aux_h__
#define __aux_h__

//corresponding to data 1
//#define nu 2
//#define n 20

// corresponding to data 2
//#define nu 3
//#define n 27

// corresponding to data 3
//#define nu 2
//#define n 10


// The basic set of data
typedef struct data data;

// The root to be constructed
typedef struct root root;

// Nodes of which the root conists
typedef struct rnode rnode;

// The points in the box
typedef struct point point;

// Initiates basic data
data * init_data(int n, int nu);

// Initiates a root
root * init_root();

// Initiates a point
point * init_point(data * I);

// Creates a new point up in the direction of s.
point * create_point(int s, point * p);

// Nifty tool
int pos(int i);

// Artin Laufer Nemethi algorithm
void ALN(int s, point * p);

// Delete a point/rnode from a linked list
void del_point(point * p);
void del_rnode(rnode * p);

// Are we still in the box if we move up i and down j?
int is_in_box(point * p, int i, int j);

// Processes root data from p up in s direction
void process_roots(point * p);

// Search till we get the ultimate owner!
rnode * ult_owner(rnode * r);

// Create a new set of rnodes going up from p int the s direction
void create_rnodes(point * p, point * u);

// Create a new dual rnode in the root R at level i
void create_drnode(int i, point * p);

// Create a new rnode in the root R at level i
void create_rnode(int i, point * p);

// Create the new roots
void finish_roots(point * p);

// Create the roots
void put_roots(point * p);

// Print the root
void print_data(data * I);

// The main procedure
void calculate_root(data * I);

// Clean up
void del_data(data * I);

// Flush out the data
void flush_root(data * I);

// Arithmetic functions
int element(int i, int s);
int intersection(int a, int b, int nup);


int Euler(root * R);
int sign(int s, int nu);
int max_chi(point * p, int s);
void get_filename(data * I);

struct root
{
	int level;
	int Euler;
	int names;
	int dnames;
	rnode * list;
	rnode * dlist;
	root * next;
};

struct rnode
{
	rnode * next;
	rnode * prev;
	int name;
	rnode * parent;
	rnode * owner;
	root * R;
};

struct data
{
	int n;
	int nu;
	int nup;
	int *m;
	int *bad;
	int *ZKbad;
	root * R;
	int min;
	rnode * unb;
	int format;
	FILE * out;
};

struct point
{
	int * coord;
	int chi;
	int min;
	int level;
	int floor;
	int ceil;
	int unbound;
	data * I;
	point * next;
	point * prev;
	point ** updown;
	rnode ** roots;
	rnode ** droots;
};
#endif
