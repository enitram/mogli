// This software has been written by Christine Solnon.
// It is distributed under the CeCILL-B FREE SOFTWARE LICENSE
// see http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for more details

#include "sublad.h"

void printDomains(Tdomain *d, int n){
	int u, i;
	for (u=0; u<n; u++){
		printf("D[%d] = ",u);
		for (i=0; i<d->nbVal[u]; i++)
			printf("%d ",d->val[d->firstVal[u]+i]);
		printf("\n");
	}
		
}

c_bool toFilterEmpty(Tdomain* D){
	// return c_true if there is no more nodes in toFilter
	return (D->nextOutToFilter < 0);
}

void resetToFilter(Tdomain *D, int size){
	// empty to filter and unmark the vertices that are marked to be filtered
	memset(D->markedToFilter,c_false,size);
	D->nextOutToFilter = -1;
}

int nextToFilter(Tdomain* D, int size){
	// precondition: toFilterEmpty = c_false
	// remove a node from toFilter (FIFO)
	// unmark this node and return it
	int u = D->toFilter[D->nextOutToFilter];
	D->markedToFilter[u] = c_false;
	if (D->nextOutToFilter == D->lastInToFilter) 
		// u was the last node in tofilter
		D->nextOutToFilter = -1;
	else if (D->nextOutToFilter == size-1)
		D->nextOutToFilter = 0;
	else D->nextOutToFilter++;
	return u;
}

void addToFilter(int u, Tdomain* D, int size){
	// if u is not marked, then add it to toFilter and mark it
	if (D->markedToFilter[u]) return;
	D->markedToFilter[u] = c_true;
	if (D->nextOutToFilter < 0){
		D->lastInToFilter = 0;
		D->nextOutToFilter = 0;
	}
	else if (D->lastInToFilter == size-1)
		D->lastInToFilter = 0;
	else D->lastInToFilter++;
	D->toFilter[D->lastInToFilter] = u;
}

c_bool isInD(int u, int v, Tdomain* D){
	// returns c_true if v belongs to D(u); c_false otherwise
	return (D->posInVal[u][v]<D->firstVal[u]+D->nbVal[u]);
}

c_bool augmentingPath(int u, Tdomain* D, int nbV){
	// return c_true if there exists an augmenting path starting from u and ending on a free vertex v
	// in the bipartite directed graph G=(U,V,E) such that U=pattern nodes, V=target nodes, and 
	// E={(u,v), v in D(u)} U {(v,u), D->globalMatchingP[u]=v}
	// update D-globalMatchingP and D->globalMatchingT consequently
	int fifo[nbV];
	int pred[nbV];
	int nextIn = 0;
	int nextOut = 0;
	int i, v, v2, u2;
	c_bool marked[nbV];
	memset(marked,c_false,nbV);
	for (i=0; i<D->nbVal[u]; i++){
		v = D->val[D->firstVal[u]+i];// v in D(u)
		if (D->globalMatchingT[v]<0){// v is free => augmenting path found
			D->globalMatchingP[u]=v;
			D->globalMatchingT[v]=u;
			return c_true;
		}
		// v is not free => add it to fifo
		pred[v] = u;
		fifo[nextIn++] = v;
		marked[v] = c_true;
	}
	while (nextOut < nextIn){
		u2 = D->globalMatchingT[fifo[nextOut++]];
		for (i=0; i<D->nbVal[u2]; i++){
			v = D->val[D->firstVal[u2]+i];// v in D(u2)
			if (D->globalMatchingT[v]<0){// v is free => augmenting path found
				while (u2 != u){// update global matching wrt path
					v2 = D->globalMatchingP[u2];
					D->globalMatchingP[u2]=v;
					D->globalMatchingT[v]=u2;
					v = v2;
					u2 = pred[v];
				}
				D->globalMatchingP[u]=v;
				D->globalMatchingT[v]=u;
				return c_true;
			}
			if (!marked[v]){// v is not free and not marked => add it to fifo
				pred[v] = u2;
				fifo[nextIn++] = v;
				marked[v] = c_true;
			}
		}
	}
	return c_false;
}

c_bool removeAllValuesButOne(int u, int v, Tdomain* D, Tgraph* Gp, Tgraph* Gt){
	// remove all values but v from D(u) and add all successors of u in toFilter
	// return c_false if an inconsistency is detected wrt to global all diff
	int j, oldPos, newPos;
	// add all successors of u in toFilter
	for (j=0; j<Gp->nbAdj[u]; j++)
		addToFilter(Gp->adj[u][j], D, Gp->nbVertices);
	// remove all values but v from D[u]
	oldPos = D->posInVal[u][v];
	newPos = D->firstVal[u];
	D->val[oldPos] = D->val[newPos];
	D->val[newPos] = v;
	D->posInVal[u][D->val[newPos]] = newPos;
	D->posInVal[u][D->val[oldPos]] = oldPos;
	D->nbVal[u] = 1;
	// update global matchings that support the global all different constraint
	if (D->globalMatchingP[u]!=v){
		D->globalMatchingT[D->globalMatchingP[u]]=-1;
		D->globalMatchingP[u]=-1;
		return augmentingPath(u,D,Gt->nbVertices);
	}
	return c_true;
}


c_bool removeValue(int u, int v, Tdomain* D, Tgraph* Gp, Tgraph* Gt){
	// remove v from D(u) and add all successors of u in toFilter
	// return c_false if an inconsistency is detected wrt global all diff
	int j;

	// add all successors of u in toFilter
	for (j=0; j<Gp->nbAdj[u]; j++)
		addToFilter(Gp->adj[u][j], D, Gp->nbVertices);
	// remove v from D[u]
	int oldPos = D->posInVal[u][v];
	D->nbVal[u]--;
	int newPos = D->firstVal[u]+D->nbVal[u];
	D->val[oldPos] = D->val[newPos];
	D->val[newPos] = v;
	D->posInVal[u][D->val[oldPos]] = oldPos;
	D->posInVal[u][D->val[newPos]] = newPos;
	// update global matchings that support the global all different constraint
	if (D->globalMatchingP[u]==v){
		D->globalMatchingP[u]=-1;
		D->globalMatchingT[v]=-1;
		return augmentingPath(u,D,Gt->nbVertices);
	}
	return c_true;
}

c_bool isCompatible(c_bool induced, int dirGp, int dirGt){
	if (dirGp == dirGt) return c_true;
	if (induced) return c_false;
	if (dirGt == 3) return c_true;
	return c_false;
}


c_bool matchVertices(int nb, int* toBeMatched, c_bool induced, Tdomain* D, Tgraph* Gp, Tgraph* Gt){
	// for each u in toBeMatched[0..nb-1], match u to D->val[D->firstVal[u]
	// and filter domains of other non matched vertices wrt FC(Edges) and FC(diff)
	// (this is not mandatory, as LAD is stronger than FC(Edges) and GAC(allDiff) 
	// is stronger than FC(diff), but this speeds up the solution process).
	// return c_false if an inconsistency is detected by FC(Edges) or FC(diff); c_true otherwise;
	int j, u, v, u2, oldNbVal;
	while (nb>0){
		u = toBeMatched[--nb];
		v = D->val[D->firstVal[u]]; 
		// match u to v
		for (u2=0; u2<Gp->nbVertices; u2++){
			if (u != u2){
				oldNbVal = D->nbVal[u2];
				if (isInD(u2,v,D) && !removeValue(u2,v,D,Gp,Gt)) return c_false;
				if (Gp->edgeDirection[u][u2] != 0){// remove from D[u2] vertices which are not adjacent to v
					j=D->firstVal[u2]; 
					while (j<D->firstVal[u2]+D->nbVal[u2]){
						if ((!Gp->isLabelled || compatibleEdgeLabels(Gp->edgeLabel[u][u2], Gt->edgeLabel[v][D->val[j]]))
							&& (isCompatible(induced, Gp->edgeDirection[u][u2], Gt->edgeDirection[v][D->val[j]]))) j++;
						else if (!removeValue(u2,D->val[j],D,Gp,Gt)) return c_false;
					}
				}
				else if (induced){// (u,u2) is not an edge => remove neighbors of v from D[u2]
					j=D->firstVal[u2]; 
					while (j<D->firstVal[u2]+D->nbVal[u2]){
						if (Gt->edgeDirection[v][D->val[j]] == 0) j++;
						else if (!removeValue(u2,D->val[j],D,Gp,Gt)) return c_false;
					}
				}
				if (D->nbVal[u2] == 0) return c_false; // D[u2] is empty
				if ((D->nbVal[u2] == 1) && (oldNbVal > 1)) toBeMatched[nb++]=u2;
			}			
		}
	}
	return c_true;
}


c_bool matchVertex(int u, c_bool induced, Tdomain* D, Tgraph* Gp, Tgraph *Gt){
	// match u to D->val[D->firstVal[u]]
	// and filter domains of other non matched vertices wrt FC(Edges) and FC(diff)
	// (this is not mandatory, as LAD is stronger than FC(Edges) and GAC(allDiff) 
	// is stronger than FC(diff), but this speeds up the solution process).
	// return c_false if an inconsistency is detected by FC(Edges) or FC(diff); c_true otherwise;
	int toBeMatched[Gp->nbVertices];
	toBeMatched[0]=u;
	return matchVertices(1,toBeMatched,induced,D,Gp,Gt);
}


int qcompare (void const *a, void const *b){
	// function used by the qsort function
	int pa = *((int*)a) - *((int*)b);
	return pa;
}

c_bool compare(int size_mu, int* mu, int size_mv, int* mv){
	// return c_true if for every element u of mu there exists
	// a different element v of mv such that u <= v; 
	// return c_false otherwise
	int i, j;
	qsort(mu, size_mu, sizeof(int), qcompare);
	qsort(mv, size_mv, sizeof(int), qcompare);
	i = size_mv-1;
	for (j=size_mu-1; j>=0; j--){
		if (mu[j]>mv[i]) return c_false;
		i--;
	}
	return c_true;
}

c_bool compatibleVertices(c_bool induced, int u, int v, Tgraph* Gp, Tgraph* Gt){
    if (Gp->isLoop[u] != Gt->isLoop[v]){
        if ((induced) || (Gp->isLoop[u]))
            return c_false;
    }
	if (Gp->isLabelled && !compatibleVertexLabels(Gp->vertexLabel[u], Gt->vertexLabel[v]))
		return c_false;
	if (!Gp->isDirected && !Gt->isDirected){ // non directed graphs
		if (Gp->nbAdj[u] > Gt->nbAdj[v]) return c_false;
		int i, mu[Gp->nbAdj[u]], mv[Gt->nbAdj[v]];
		for (i=0; i<Gp->nbAdj[u]; i++) mu[i]=Gp->nbAdj[Gp->adj[u][i]];
		for (i=0; i<Gt->nbAdj[v]; i++) mv[i]=Gt->nbAdj[Gt->adj[v][i]];
		return compare(Gp->nbAdj[u], mu, Gt->nbAdj[v], mv);
	}
	// directed graphs
	if (induced){
		if (Gp->nbPred[u] > Gt->nbPred[v]) return c_false;
		if (Gp->nbSucc[u] > Gt->nbSucc[v]) return c_false;
		if (Gp->nbAdj[u] - Gp->nbPred[u] - Gp->nbSucc[u] > Gt->nbAdj[v] - Gt->nbPred[v] - Gt->nbSucc[v]) return c_false;
		return c_true;
	}
	if (Gp->nbAdj[u] > Gt->nbAdj[v]) return c_false;
	if (Gp->nbPred[u] > Gt->nbAdj[v] - Gt->nbSucc[v]) return c_false;
	if (Gp->nbSucc[u] > Gt->nbAdj[v] - Gt->nbPred[v]) return c_false;
	if (Gp->nbAdj[u] - Gp->nbPred[u] - Gp->nbSucc[u] > Gt->nbAdj[v] - Gt->nbPred[v] - Gt->nbSucc[v]) return c_false;
	return c_true;
}

c_bool initDomains(c_bool induced, Tdomain* D, Tgraph* Gp, Tgraph* Gt){
	// for every pattern node u, initialize D(u) with every vertex v 
	// such that for every neighbor u' of u there exists a different 
	// neighbor v' of v such that degree(u) <= degree(v)
	// if initialDomains, then filter initial domains wrt compatibilities given in file
	// return c_false if a domain is empty and c_true otherwise
	int val[Gp->nbVertices*Gt->nbVertices];
	int matchingSize, u, v, i, valSize;
	matchingSize = 0;
	valSize = 0;
	for (u=0; u<Gp->nbVertices; u++){
		D->markedToFilter[u] = c_true;
		D->toFilter[u] = u;
		D->nbVal[u] = 0;
		D->posInVal[u] = (int*)malloc(sizeof(int)*Gt->nbVertices);
		D->firstMatch[u] = (int*)malloc(sizeof(int)*Gt->nbVertices);
		D->firstVal[u] = valSize;
		for (v=0; v<Gt->nbVertices; v++){
			if (!compatibleVertices(induced, u, v, Gp, Gt)) // v not in D(u)
				D->posInVal[u][v] = D->firstVal[u]+Gt->nbVertices;
			else { // v in D[u]
				D->firstMatch[u][v] = matchingSize;
				matchingSize += Gp->nbAdj[u];
				val[valSize] = v;
				D->nbVal[u]++;
				D->posInVal[u][v] = valSize++;
			}
		}
		if (D->nbVal[u]==0) return 0; // empty domain
	}
	D->val = (int*)malloc(sizeof(int)*valSize);
	for (i=0; i<valSize; i++) D->val[i] = val[i];
	D->matching = (int*)malloc(sizeof(int)*matchingSize);
	memset(D->matching,-1,sizeof(int)*matchingSize);
	D->nextOutToFilter = 0;
	D->lastInToFilter = Gp->nbVertices-1;
	return 1;
}

