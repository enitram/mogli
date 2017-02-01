// This software has been written by Christine Solnon.
// It is distributed under the CeCILL-B FREE SOFTWARE LICENSE
// see http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for more details

#include "sublad.h"

// Global variables
int nbNodes = 1;      // number of nodes in the search tree
int nbFail = 0;       // number of failed nodes in the search tree
int nbSol = 0;        // number of solutions found
struct rusage ru;     // reusable structure to get CPU time usage

c_bool filter(c_bool induced, Tdomain* D, Tgraph* Gp, Tgraph* Gt){
	// filter domains of all vertices in D->toFilter wrt LAD and ensure GAC(allDiff)
	// return false if some domain becomes empty; true otherwise
	int u, v, i, oldNbVal;
	while (!toFilterEmpty(D)){
		while (!toFilterEmpty(D)){
			u=nextToFilter(D,Gp->nbVertices); 
			oldNbVal = D->nbVal[u];
			i=D->firstVal[u]; 
			while (i<D->firstVal[u]+D->nbVal[u]){
				// for every target node v in D(u), check if G_(u,v) has a covering matching
				v=D->val[i]; 
				if (checkLAD(induced,u,v,D,Gp,Gt)) i++;
				else if (!removeValue(u,v,D,Gp,Gt)) return c_false;
			}
			if ((D->nbVal[u]==1) && (oldNbVal>1) && (!matchVertex(u,induced,D,Gp,Gt))) return c_false;
			if (D->nbVal[u]==0) return c_false;
		}
		if (!ensureGACallDiff(induced,Gp,Gt,D)) return c_false;
	}
	return c_true;
}



c_bool solve(int timeLimit, c_bool firstSol, c_bool induced, int verbose, Tdomain* D, Tgraph* Gp, Tgraph* Gt, int *map){
	// if firstSol then search for the first solution; otherwise search for all solutions
	// if induced then search for induced subgraphs; otherwise search for partial subgraphs
	// return false if CPU time limit exceeded before the search is completed
	// return true otherwise
	
	int u, v, minDom, i; 
	int nbVal[Gp->nbVertices];
	int globalMatching[Gp->nbVertices];
	
	nbNodes++;

	getrusage(RUSAGE_SELF, &ru);
	if (ru.ru_utime.tv_sec >= timeLimit)	
		// CPU time limit exceeded
		return c_false;
	
	if (!filter(induced,D,Gp,Gt)){ 
		// filtering has detected an inconsistency
		if (verbose == 2) printf("Filtering has detected an inconsistency\n");
		nbFail++;
		resetToFilter(D,Gp->nbVertices);
		return c_true;
	}	
	
	// The current node of the search tree is consistent wrt to LAD and GAC(allDiff)
	// Save domain sizes and global all different matching
	// and search for the non matched vertex minDom with smallest domain
	minDom=-1;
	for (u=0; u<Gp->nbVertices; u++){
		nbVal[u]=D->nbVal[u];
		if ((nbVal[u]>1) && ((minDom<0) || (nbVal[u]<nbVal[minDom]))) minDom=u;
		globalMatching[u] = D->globalMatchingP[u];
	}
	
	if (minDom==-1){ 
		// All vertices are matched => Solution found
		nbSol++;
    for (u=0; u<Gp->nbVertices; u++) {
      map[u] = D->val[D->firstVal[u]];
    }
		if (verbose >= 1){
			printf("Solution %d: ",nbSol);
			for (u=0; u<Gp->nbVertices; u++) {
        printf("%d=%d ",u,D->val[D->firstVal[u]]);
      }
			printf("\n");
		}
		resetToFilter(D,Gp->nbVertices);
		return c_true;
	}
	
	// save the domain of minDom to iterate on its values
	int val[D->nbVal[minDom]];
	for (i=0; i<D->nbVal[minDom]; i++) val[i]=D->val[D->firstVal[minDom]+i];
	
	// branch on minDom=v, for every target node v in D(u)
	for(i=0; ((i<nbVal[minDom]) && ((firstSol==0)||(nbSol==0))); i++){
		v = val[i]; 
		if (verbose == 2) printf("Branch on %d=%d\n",minDom,v);
		if ((!removeAllValuesButOne(minDom,v,D,Gp,Gt)) || (!matchVertex(minDom,induced,D,Gp,Gt))){
			if (verbose == 2) printf("Inconsistency detected while matching %d to %d\n",minDom,v);
			nbFail++; 
			nbNodes++;
			resetToFilter(D,Gp->nbVertices);
		} 
		else if (!solve(timeLimit,firstSol,induced,verbose,D,Gp,Gt,map))
			// CPU time exceeded
			return c_false;
		// restore domain sizes and global all different matching
		if (verbose == 2) printf("End of branch %d=%d\n",minDom,v);
		memset(D->globalMatchingT,-1,sizeof(int)*Gt->nbVertices);
		for (u=0; u<Gp->nbVertices; u++){
			D->nbVal[u] = nbVal[u];
			D->globalMatchingP[u] = globalMatching[u];
			D->globalMatchingT[globalMatching[u]] = u;
		}
	}
	return c_true;
}

int printStats(c_bool timeout, int verbose){
  if (verbose > 0) {
    // print statistics line and return exit status depending on timeout
    getrusage(RUSAGE_SELF, &ru);
    if (timeout)
      printf("CPU time exceeded");
    else
      printf("Run completed");
    printf(": %d solutions; %d fail nodes; %d nodes; %d.%06d seconds\n",
           nbSol, nbFail, nbNodes,
           (int) ru.ru_utime.tv_sec, (int) ru.ru_utime.tv_usec);
  }
  return timeout;
}

int sub_iso(Tgraph *Gp, Tgraph *Gt, c_bool *iso, int *map, int verbose, int timeLimit, c_bool induced, c_bool firstSol) {
	// Parameters
//	int timeLimit=60;      // Default: CPU time limit set to 60 seconds
//	int verbose = 0;       // Default: non verbose execution
//	c_bool induced = c_true;  // Default: search for partial subgraph
//	c_bool firstSol = c_true; // Default: search for all solutions
  *iso = c_false;

	nbNodes = 1;      // number of nodes in the search tree
	nbFail = 0;       // number of failed nodes in the search tree
	nbSol = 0;        // number of solutions found

  matchedWithV = (int*)malloc(Gt->nbVertices*sizeof(int));
  nbPred = (int*)malloc(Gt->nbVertices*sizeof(int));
  pred = (int**)malloc(Gt->nbVertices*sizeof(int*));
  succ  = (int**)malloc(Gt->nbVertices*sizeof(int*));
  for (int i=0; i<Gt->nbVertices; i++){
      pred[i] = (int*)malloc(Gt->nbVertices*sizeof(int));
      succ[i] = (int*)malloc(Gt->nbVertices*sizeof(int));
  }
  nbSucc  = (int*)malloc(Gt->nbVertices*sizeof(int));
  listV  = (int*)malloc(Gt->nbVertices*sizeof(int));
  listDV  = (int*)malloc(Gt->nbVertices*sizeof(int));
  listU  = (int*)malloc(Gp->nbVertices*sizeof(int));
  listDU  = (int*)malloc(Gp->nbVertices*sizeof(int));

	// Initialize domains
	Tdomain* D = (Tdomain*)malloc(sizeof(Tdomain));
	D->globalMatchingP = (int*)malloc(sizeof(int)*Gp->nbVertices);
	memset(D->globalMatchingP,-1,sizeof(int)*Gp->nbVertices);
	D->globalMatchingT = (int*)malloc(sizeof(int)*Gt->nbVertices);
	memset(D->globalMatchingT,-1,sizeof(int)*Gt->nbVertices);
	D->nbVal = (int*)malloc(sizeof(int)*Gp->nbVertices);
	D->firstVal = (int*)malloc(sizeof(int)*Gp->nbVertices);
	D->posInVal = (int**)malloc(sizeof(int*)*Gp->nbVertices);
	D->firstMatch = (int**)malloc(sizeof(int*)*Gp->nbVertices);
	D->markedToFilter = (c_bool*)calloc(Gp->nbVertices,sizeof(c_bool));
	D->toFilter = (int*)malloc(sizeof(int)*Gp->nbVertices);

	if (!initDomains(induced, D, Gp, Gt)) return printStats(c_false, verbose);
	if (verbose >= 2) printDomains(D, Gp->nbVertices);

	// Check the global all different constraint                                                                                                    
	if ((!updateMatching(Gp->nbVertices,Gt->nbVertices,D->nbVal,D->firstVal,D->val,D->globalMatchingP)) ||
		(!ensureGACallDiff(induced,Gp,Gt,D))){
		nbFail++;
		return printStats(c_false, verbose);
	}

	// Math all vertices with singleton domains                                                                                                     
	int u;
	int nbToMatch = 0;
	int toMatch[Gp->nbVertices];
	for (u=0; u<Gp->nbVertices; u++){
		D->globalMatchingT[D->globalMatchingP[u]] = u;
		if (D->nbVal[u] == 1)
			toMatch[nbToMatch++] = u;
	}
	if (!matchVertices(nbToMatch,toMatch,induced,D,Gp,Gt)){
		nbFail++;
		return printStats(c_false, verbose);
	}

	// Solve
  c_bool timeout = !solve(timeLimit,firstSol, induced, verbose, D, Gp, Gt, map);
  if (nbSol > 0) {
    *iso = c_true;
  }
	int ret = printStats(timeout, verbose);

	free(matchedWithV);
	free(nbPred);
	for (int i=0; i<Gt->nbVertices; i++){
		free(pred[i]);
		free(succ[i]);
	}
	free(pred);
	free(succ);
	free(nbSucc);
	free(listV);
	free(listDV);
	free(listU);
	free(listDU);

	free(D->globalMatchingP);
	free(D->globalMatchingT);
	free(D->nbVal);
	free(D->firstVal);
	free(D->posInVal);
	free(D->firstMatch);
	free(D->markedToFilter);
	free(D->toFilter);
	free(D);

	return ret;
}

