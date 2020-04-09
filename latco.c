#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#include "latco.h"


// The basic set of data to be analized

data * init_data(int n, int nu)
{
	data * I = malloc(sizeof(data));
	I->n = n;
	I->nu = nu;
	I->nup = 1<<nu;
	I->bad = (int*) malloc(nu*sizeof(int));
	I->ZKbad = (int*) malloc(nu*sizeof(int));
	I->m = (int*) malloc(n*n*sizeof(int));
	I->min = 0;
	I->R = NULL;
	I->unb = (rnode*) malloc(sizeof(rnode));
	I->unb->owner = NULL;
	I->format = 0;
	I->out = NULL;
	return I;
}

root * init_root()
{
	root * R = malloc(sizeof(root));
	R->level = 0;
	R->Euler = -1;
	R->names = 0;
	R->next = NULL;
	R->list = NULL;
    return R;

//	R->dlist = NULL;
//	R->dlist->next = NULL;
//	R->dlist->prev = NULL;
//	R->dlist->name = R->dnames++;
//	R->dlist->parent = NULL;
//	R->dlist->owner = NULL;
//	R->dlist->R = R;
}

point * init_point(data * I)
{
	int i,j;
	point * p = malloc(sizeof(point));
	p->coord = malloc(I->n*sizeof(int));
	for(i=0; i<I->n; i++)
		p->coord[i] = 0;
	p->chi = 0;
//	p->dchi = 0;
	p->min = 0;
	p->level = 0;
	p->I = I;
	p->next = NULL;
	p->prev = NULL;
	p->floor = I->nup - 1;
	p->ceil = 0;
	p->unbound = 0;
	p->updown = malloc(I->nup*I->nup*sizeof(point*));
	for(i=0; i<I->nup; i++)
	{
		for(j=0; j<I->nup; j++)
		{
			if(i!=j)
				p->updown[I->nup*i + j] = NULL;
			else
				p->updown[I->nup*i + i] = p;
		}
	}
	p->roots = NULL;
	p->droots = NULL;
	return p;
}

point * create_point(int s, point * p)
{
	int i,j;
	int nup = p->I->nup;
	point * new = malloc(sizeof(point));
	new->roots = NULL;
	new->droots = NULL;
	new->coord = malloc(p->I->n*sizeof(int));
	new->I = p->I;

	int sm = 0;
	while( (s & (1<<sm)) == 0)
		sm++;
	point * below = p->updown[nup*(s - (1<<sm))];
	new->unbound = p->unbound;

	new->level = 0;

	for(i=0; i<new->I->n; i++)
	{
		new->coord[i] = below->coord[i];
	}

	new->chi = below->chi;
//	new->dchi = 0;

	ALN(sm, new);
	new->min = below->min + pos(below->chi - new->chi);

	new->ceil = 0;
	new->floor = 0;
	for(i=0; i<p->I->nu; i++)
	{
		if(new->coord[new->I->bad[i]] == 0)
			new->floor += (1<<i);
		if(new->coord[new->I->bad[i]] == new->I->ZKbad[i])
			new->ceil += (1<<i);
		new->level += new->coord[new->I->bad[i]];
	}

	// This configures updown data pointing at or from the new point
	
	new->updown = malloc(p->I->nup*p->I->nup*sizeof(point*));
	int ip, jp, c;

	for(i=0; i<nup; i++)
	{
		new->updown[nup*i + i] = new;
		for(j=0; j<nup; j++)
		{
			ip = i - (i&j);
			jp = j - (i&j);
			if((ip&new->floor) || (jp&new->ceil) || (jp==nup-1))
			{
				new->updown[nup*j + i] = NULL;
				continue;
			}
			if(i == j)
			{
				new->updown[nup*j + i] = new;
				continue;
			}

			c = 0;
			while((1<<c) & (jp | new->floor))
				c++;
			if(c == new->I->nu)
			{
				new->updown[nup*j + i] = NULL;
				continue;
			}

			new->updown[nup*j + i]
				= p->updown[nup*s           + (1<<c)]
				   ->updown[nup*(jp+(1<<c)) +     ip];
			if(new->updown[nup*j + i] != NULL)
				new->updown[nup*j + i]->updown[nup*i+j] = new;
		}
	}



//	printf("created  %d %d -- %p -- level = %d\n",
//		new->coord[new->I->bad[0]],
//	    new->coord[new->I->bad[1]], new, new->level);

	return new;
}

int is_in_box(point * p, int i, int j)
{
	if(p->coord[p->I->bad[i]] < p->I->ZKbad[i]
		&& 0 < p->coord[p->I->bad[j]])
		return 1;
	return 0;
}

int pos(int i)
{
	if(i > 0)
		return i;
	else
		return 0;
}

void ALN(int s, point * p)
{
	int i, j;
	int d;
	int Bad[p->I->n];
	int dummy1 = 1;
	int dummy2 = 1;

	for(i=0; i<p->I->n; i++)
		Bad[i] = 0;

	for(i=0; i<p->I->nu; i++)
		Bad[p->I->bad[i]] = 1;



	// The first step in the sequence (the one that provides possible
	// contribution to the geometric genus!)
	d = 0;
	for(j=0; j<p->I->n; j++)
		d += p->coord[j] * p->I->m[p->I->bad[s]*p->I->n+j];
	p->chi += 1-d;
	p->coord[p->I->bad[s]]++;
	p->level++;
	
	// The rest of them
	while(dummy1)
	{
		dummy1 = 0;
		dummy2 = 1;
		for(i=0; i<p->I->n && dummy2; i++)
		{
			if(!Bad[i])
			{
				d = 0;
				for(j=0; j<p->I->n; j++)
					d += p->coord[j] * p->I->m[p->I->n*i+j];
				if(d > 0)
				{
					p->chi += 1 - d;
					p->coord[i]++;
					dummy2 = 0;
					dummy1 = 1;
				}
			}
		}
	}
}

void del_point(point * p)
{
	free(p->coord);
	free(p->updown);
	free(p->roots);
	free(p->droots);
	if(p->next != NULL)
		p->next->prev = p->prev;
	if(p->prev != NULL)
		p->prev->next = p->next;
	free(p);
}

void process_roots(point * p)
{
	int i,j;
	if(p->chi<=0)
		p->roots = (rnode**) malloc((-p->chi+1)*sizeof(rnode*));
	else
		p->roots = NULL;
	point * d;
	rnode * dow;
//    root * R = p->I->R;
	for(i=0; i>=p->chi; i--)
	{
		p->roots[-i] = NULL;
		for(j=0; j<p->I->nu; j++)
		{
			if(p->floor & (1<<j))
				continue;

			d = p->updown[1<<j];
			if(d->chi <= i)
			{
				dow = ult_owner(d->roots[-i]);
				if(p->roots[-i] == NULL)
				{
					p->roots[-i] = dow;
				}
				else
				{
					if(dow != p->roots[-i])
						dow->owner = p->roots[-i];
				}
			}
		}
		if(p->roots[-i] == NULL)
			create_rnode(i,p);
	}



	for(j=1; j<p->I->nup; j++)
		if(!(j & p->floor) && p->updown[j]->unbound > p->unbound)
			p->unbound = p->updown[j]->unbound;

	if(p->unbound > p->chi-1 || p->floor || p->ceil)
	{
		p->unbound = p->chi-1;
	}

	int nrchi;

	if(p->chi - p->unbound-1 > 0)
		p->droots = (rnode**) malloc((p->chi - p->unbound-1) * sizeof(rnode*));

	if(p->chi-1 != p->unbound)
		printf("hello %d %d  %d %d\n", p->coord[p->I->bad[0]], p->coord[p->I->bad[1]], p->chi, p->unbound);
	for(i=p->chi-1; i>p->unbound; i--)
	{
		p->droots[-i +p->chi-1] = NULL;
		for(j=1; j<p->I->nup; j++)
		{
			d = p->updown[j];
			nrchi = -pos(-d->chi);
			dow = ult_owner(d->droots[-i +nrchi-1]);
			if(nrchi > i && i > d->unbound)
			{
				if(p->droots[-i +p->chi-1] == NULL)
					p->droots[-i +p->chi-1] = dow;
				else if(p->droots[-i +p->chi-1] == p->I->unb
				                         && dow != p->I->unb)
				 	dow->owner = p->I->unb;
				else if(p->droots[-i +p->chi-1] != dow)
					ult_owner(p->droots[-i +p->chi-1])->owner = dow;
			}
		}
		if(p->droots[-i +p->chi-1] == NULL)
		{
			create_drnode(i,p);
		}
	}

	for(j=1; j<p->I->nup; j++)
	{
		if(!(j & p->floor))
		{
			d = p->updown[j];
			nrchi = -pos(-d->chi);
			if(p->unbound < nrchi-1)
				i = p->unbound;
			else
				i = nrchi-1;
			for(; i>d->unbound; i--)
			{
				dow = ult_owner(d->droots[-i +nrchi-1]);
				if(dow != p->I->unb)
					dow->owner = p->I->unb;
			}
		}
	}
	
//	p->droots = (rnode*) malloc(
/*	if(p->dchi > 0)
		return;
	p->droots = (rnode**) malloc((-p->chi + p->dchi)*sizeof(rnode*));
	for(i=p->dchi-1; i>=p->chi; i--)
	{
		p->droots[-i +p->dchi-1] = NULL;
		for(j=0; j<p->I->nu; j++)
		{
			if(p->floor & (1<<j))
			{
				p->droots[-i +p->dchi-1] = p->roots[-i]->R->dlist;
				continue;
			}
			d = p->updown[1<<j];
			if(d->dchi > 0)
			{
				p->droots[-i +p->dchi-1] = p->roots[-i]->R->dlist;
				continue;
			}
			if(d->dchi > i && d->chi <= i)
			{
				dow = ult_owner(d->droots[-i +d->dchi-1]);
				if(p->droots[-i +p->dchi-1] == NULL)
					p->droots[-i +p->dchi-1] = dow;
				else if(dow == p->roots[-i]->R->dlist
				        && p->droots[-i +p->dchi-1] != p->roots[-i]->R->dlist)
					p->droots[-i +p->dchi-1]->owner = dow;
				else if(p->droots[-i +p->dchi-1] != dow)
					dow->owner = p->droots[-i +p->dchi-1];
			}
		}
		if(p->droots[-i +p->dchi-1] == NULL)
			create_drnode(i,p);
	}*/
}

void del_rnode(rnode * rn)
{
	if(rn->next != NULL)
		rn->next->prev = rn->prev;
	if(rn->prev != NULL)
		rn->prev->next = rn->next;
	free(rn);
}

rnode * ult_owner(rnode * r)
{
	if(r == NULL)
	{
		printf("Can not compute ult_owner(NULL).\n");
		return NULL;
	}
	if(r->owner == r)
	{
		printf("rnode cannot be its own owner!\n");
		return NULL;
	}
	while(r->owner != NULL)
	{
		r = r->owner;
	}
	return r;
}

void create_drnode(int i, point * p)
{
	printf("Creating rnode!\n");
	root * R = p->I->R;
	while(R->level != i)
		R = R->next;
	p->droots[-i +p->chi-1] = (rnode*) malloc(sizeof(rnode));

	p->droots[-i +p->chi-1]->next = R->dlist;
	if(p->droots[-i +p->chi-1]->next != NULL)
		p->droots[-i +p->chi-1]->next->prev = p->droots[-i +p->chi-1];
	p->droots[-i +p->chi-1]->prev = NULL;

	p->droots[-i +p->chi-1]->owner = NULL;

	p->droots[-i +p->chi-1]->name = R->dnames++;
}

void create_rnode(int i, point * p)
{
	if(i < p->chi || i>0)
	{
		printf("This is no good\n");
		return;
	}
	root * R = p->I->R;
	rnode * new = malloc(sizeof(rnode));
	while(R->level != i)
	{
		if(R->next == NULL)
		{
//			printf("New level: %d.\n", R->level - 1);
			R->next = malloc(sizeof(root));
			R->next->level = R->level - 1;
			R->next->Euler = -1;
			R->next->names = 0;
			R->next->dnames = 0;
			R->next->list = NULL;
			R->next->dlist = NULL;
			R->next->next = NULL;
		}
		R = R->next;
	}
	new->next = R->list;
	new->prev = NULL;
	if(R->list != NULL)
		R->list->prev = new;
	R->list = new;
	new->R = R;
	new->name = R->names++;

	if(i!=0)
		new->parent = ult_owner(p->roots[-i-1]);
	else
		new->parent = NULL;
	
	new->owner = NULL;
	p->roots[-i] = new;
}

void print_data(data * I)
{
	if(I->R == NULL)
	{
		printf("There is no root!\n");
		return;
	}
	root * run = I->R;
	rnode * rnrun = NULL;
	int rkn = 0;
	int lev = 1;
	while(run != NULL)
	{
		printf("\nLevel %d:\n", run->level);
		printf("Reduced Euler characteristic is %d.\n", run->Euler);
		rnrun = run->list;
		rkn--;
		while(rnrun != NULL)
		{
			if(rnrun->owner == NULL)
			{
				printf("Component named %d", rnrun->name);
				if(rnrun->parent == NULL)
					printf(" has the unique parent.\n");
				else
					printf(" has parent %d.\n", ult_owner(rnrun->parent)->name);
				rkn++;
			}
			rnrun = rnrun->next;
		}

/*		rnrun = run->dlist;
		printf("\n");
		while(rnrun != NULL)
		{
			if(rnrun->owner == NULL)
			{
				printf("Complementary component named %d", rnrun->name);
				if(rnrun->parent == NULL)
					printf(" has the unique parent.\n");
				else
					printf(" has parent %d.\n", ult_owner(rnrun->parent)->name);
			}
			rnrun = rnrun->next;
		} */



		run = run->next;
		lev--;
	}

	printf("\n");
	printf("eu(H^*):_____________________ %d\n", Euler(I->R) - lev);
	printf("Minimal value of chi:________%d\n", lev);
	printf("rk(H^0_red):_________________ %d\n", rkn);
	printf("eu(H^0):_____________________ %d\n", rkn - lev);
	printf("Minimal path cohomology:_____ %d\n",
		I->min);
	printf("WARNING: The minimal path cohomology presented here is only an upper bound.\n");
	printf("Counterexamples are not known, but a reduction theorem in the style of László\n");
	printf("has not been proved for this invariant.\n");
}


void calculate_root(data * I)
{
	root * R = init_root();
	I->R = R;
	point * p = init_point(I);

	int s,i,dummy,counter,neg;

	point * step;
	point * bot = p;
	point * top = p;
	point * run;
	point * bot1;

    for(i=0; i<I->nu; i++)
    {
        if(I->ZKbad[i] < 0)
            return;
        else if(I->ZKbad[i] == 0)
            p->ceil += 1<<i;
    }



	if(I->format)
		I->out = fopen("output", "w");

	if(I->out == NULL)
		I->format = 0;
	
	int total = 1;
	for(i=0; i<I->nu; i++)
		total *= I->ZKbad[i];

	counter = 0;
	neg = 0;
	int sg, m;
	root * rd;

	printf("\nNumber of points is about %d million.\n", total/1000000);
	while(p != NULL)
	{
		counter++;
		if(p->chi<1)
			neg++;

		if(counter % 1000000 == 0)
			printf("%d / %d\n", counter/1000000, total/1000000);

		// The following constructs all relevant points in the cube with
		// smallest vertex p
		for(s=1; s<I->nup; s++)
		{
			if(!(p->ceil & s))
			{
				step = p->updown[p->I->nup*s];
				if(step == NULL)
				{
					step = create_point(s, p);
					
					run = top;
					while(run->level > step->level)
						run = run->prev;
					
					step->next = run->next;
					step->prev = run;
					if(step->next != NULL)
						step->next->prev = step;
					if(step->prev != NULL)
						step->prev->next = step;
					if(run == top)
						top = step;
				}
			}
		}
		process_roots(p);
		
		for(i=0; i<I->nu; i++)
		{
			if(!(p->floor&(1<<i)))
			{
				dummy = p->updown[1<<i]->min
				         + pos(p->updown[1<<i]->chi - p->chi);
				if(p->min > dummy)
				p->min = dummy;
			}
		}

		for(s=0; s<I->nup; s++)
		{
			if(!(p->ceil & s))
			{
				rd = I->R;
				sg = sign(s, I->nu);
				m = max_chi(p, s);
				while(rd != NULL && rd->level >= m)
				{
					rd->Euler += sg;
					rd = rd->next;
				}
			}
		}
		
		while(bot->level < p->level - p->I->nu)
		{
			bot1 = bot->next;
			del_point(bot);
			bot = bot1;
		}

		I->min = p->min;

		if(p->I->format)
		{
			fprintf(p->I->out, "m[%d", p->coord[p->I->bad[0]] + 1);
			for(i=1; i<p->I->nu; i++)
				fprintf(p->I->out, ",%d", p->coord[p->I->bad[i]] + 1);
			fprintf(p->I->out, "] = %d\n", p->chi);
		}

		p = p->next;
	}
	if(I->format)
		fclose(I->out);

	while(bot != NULL)
	{
		bot1 = bot->next;
		del_point(bot);
		bot = bot1;
	}
}

int Euler(root * R)
{
	int r = 0;
	while(R != NULL)
	{
		r += R->Euler;
		R = R->next;
	}
	return r;
}

int sign(int s, int nu)
{
	int i;
	int r = 0;
	for(i=0; i<nu; i++)
		if(s&(1<<i))
			r++;
	if(r/2*2 == r)
		return  1;
	else
		return -1;
}

int max_chi(point * p, int s)
{
	int r = p->chi;
	int i;
	for(i=0; i<p->I->nup; i++)
		if((i&s) == i && p->updown[p->I->nup*i] != NULL
		              && p->updown[p->I->nup*i]->chi > r)
			r = p->updown[p->I->nup*i]->chi;
	return r;
}

void del_data(data * I)
{
	free(I->m);
	free(I->bad);
	free(I->ZKbad);
	flush_root(I);
	free(I->unb);
	free(I);
}

void flush_root(data * I)
{
	root * R = I->R;
	root * R1;
	rnode * rn;
	rnode * rn1;
	while(R != NULL)
	{
		rn = R->list;
		while(rn != NULL)
		{
			rn1 = rn->next;
			del_rnode(rn);
			rn = rn1;
		}
		R1 = R->next;
		free(R);
		R = R1;
	}
	I->R = NULL;
	I->min = 0;
}

void get_filename(data * I)
{
	/*
	 * Since I don't know how to handle filenames, this function doesnt't
	 * do anything right now.
	 */
}
