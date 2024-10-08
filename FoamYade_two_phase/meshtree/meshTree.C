// Deepak Kunhappan (deepak.Kunhappan@3sr-grenoble.fr, deepak.kn1990@gmail.com)
// Part of YADE - OpenFOAM coupling, Build a  k-d tree from the mesh, and
// search for cell centers within a given distance range 
// Inspired from the C version at https://rosettacode.org/wiki/K-d_tree#C
// Header file  :meshTree.H 

#include "meshTree.H"

using namespace Foam;

void Foam::meshTree::build_tree() 
{
  std::vector<meshpt>  ptlist; 
  // label cellI;  
  forAll(meshm, cellI)
    ptlist.push_back(meshpt(&meshm[cellI].x(), &meshm[cellI].y(), &meshm[cellI].z(), cellI));
  numlevels =0; 
  root = recursive_build_tree(ptlist, 0);
}

kdNode* Foam::meshTree::recursive_build_tree(std::vector<meshpt>& ptlist, int depth) 
{

  if (! ptlist.size()){return NULL; } 
  numlevels = numlevels+1; 
  int axis = depth%ndim; 
 
  get_median(ptlist, axis);
  vec_sz md = ptlist.size()/2;  
  kdNode* node = new kdNode(ptlist[md]);
  
  std::vector<meshpt> pv1 = std::vector<meshpt>(ptlist.begin(), ptlist.begin()+md); 
  std::vector<meshpt> pv2 = std::vector<meshpt>(ptlist.begin()+md+1, ptlist.end()); 

  node->left = recursive_build_tree(pv1, depth+1);
  node->right = recursive_build_tree(pv2, depth+1);  

  return node;  
}

void Foam::meshTree::traversTree()
{
  _traversTree(root); 

}


void Foam::meshTree::get_median(std::vector<meshpt>& ptlist, const int& axis) 
{
//  std::sort(ptlist.begin(), ptlist.end(), cmpvec(axis)); too slow 
 vec_sz md = ptlist.size()/2;   
  std::nth_element(ptlist.begin(), ptlist.begin()+md, ptlist.end(), cmpvec(axis)); 
}


double Foam::meshTree::distance(const meshpt& p2, const meshpt& p1)
{
  double dist = 0.0; 
  for (unsigned int i=0; i!= p1.pt.size();++i)
  {
    double ds = *(p1.pt[i])-*(p2.pt[i]); 
    dist += ds*ds; 
  }
  return dist; 

}

int Foam::meshTree::nearestCell(const vector& p)
{
  
//  kdNode* bnode;  
  meshpt v(& p.x(), &p.y(), &p.z(), -1);
     
  double dist = distance(root->p,  v); 
  kdNode* bnode =recursive_nearest_cell(root, v, root,dist,0); 
  dist = distance(bnode->p, v); 
  dist = Foam::sqrt(dist); 
  
  return bnode->p.id; 
}

//
kdNode* Foam::meshTree::recursive_nearest_cell(kdNode* node, const meshpt& v, kdNode* best, double& best_dist,int depth)
{

  if (node ==NULL) return NULL;  
  kdNode* best1=best;  
  double dist_l = best_dist; 
  
  double distsq = distance(node->p, v); 

   if (distsq < best_dist){ 
      dist_l = distsq;
      best1 = node; 
   }

   int axis = depth%ndim;
   double df = *(node ->p.pt[axis])-*(v.pt[axis]);
   double df2 = df*df; 
   
//    std::cout << "best dist = " << dist_l  << std::endl;
//    std::cout << "df = " << df2  << std::endl;

   kdNode* next; 
   kdNode* other; 

   if (df > 0.0){
     next = node->left; 
     other =node->right;
   } else { 
     next = node-> right; 
     other = node-> left;
   }
  
  depth = depth+1;
  kdNode* nextN = recursive_nearest_cell(next, v, best1, dist_l, depth); 
  if (nextN != NULL){
    distsq =distance(nextN->p, v); 
    if (distsq < dist_l){ 
      dist_l = distsq;
      best1 = nextN;  
  }
}

  if(df2 < dist_l ){  
    kdNode* nextN = recursive_nearest_cell(other, v,best1, dist_l, depth); 
  if (nextN != NULL){
    distsq =distance(nextN->p, v); 
    if (distsq < dist_l){ 
      dist_l = distsq;
      best1 = nextN;  
    }
  }

}
return best1;
}


void Foam::meshTree::_traversTree(kdNode* node) 
{
 
  if(node==NULL){return; }
  _traversTree(node->left); 
  _traversTree(node ->right);   

}


std::vector<int> Foam::meshTree::nnearestCellsRange(const vector& v, const double& range, const bool& np)

{

  meshpt px(&v.x(), &v.y(), &v.z(), -1); 
  const unsigned int maxelem = 12; 
  pqueue pq(maxelem);  // = new pqueue(maxelem);
  pq.maxdist = (range*range )+(0.25*range*range);  
  double dist = distance(root->p, px); 
  kdNode* bestnode = nnearest(root,px, root, dist,0, pq);

  std::vector<int> ids; //if (!pq -> container.size()) ids.push_back(bestnode -> p.id); 

  if ( !np || (pq.container.size()==0)) {ids.push_back(bestnode -> p.id);}
  else {
  for (unsigned int i=0; i!= pq.container.size(); ++i)
  {
    int id = pq.container[i].first->p.id; 
    ids.push_back(id); 
    
  } }

//  for (unsigned int i=0; i != pq -> container.size(); ++ i) {
//      delete pq -> container[i].first; 
//  }

  pq.container.clear(); 
// delete pq;
//   delete bestnode; 
  return ids; 

}


kdNode* Foam::meshTree::nnearest(kdNode* node, const meshpt& v, kdNode* best, double& best_dist, int depth, pqueue& pq) {

 
  if (node ==NULL) return NULL;  

  kdNode* best1=best;  
  double dist_l = best_dist; 
  
  double distsq = distance(node->p, v); 

   if (distsq < best_dist){ 
      dist_l = distsq;
      best1 = node;
      if (dist_l < pq.maxdist) 
        pq.push_node(std::make_pair(best1, dist_l)); 
   }

   int axis = depth%ndim;
   double df = *(node ->p.pt[axis])-*(v.pt[axis]);
   double df2 = df*df; 
   
   kdNode* next; 
   kdNode* other; 

   if (df > 0.0){
     next = node->left; 
     other =node->right;
   } else { 
     next = node-> right; 
     other = node-> left;
   }
  depth = depth+1;
  kdNode* nextN  = nnearest(next, v, best1, dist_l, depth, pq); 
  if (nextN != NULL){
    distsq =distance(nextN->p, v); 
    if (distsq < dist_l){ 
      dist_l = distsq;
      best1 = nextN;
      if (dist_l < pq.maxdist) 
        pq.push_node(std::make_pair(best1, dist_l));      
  }
}

  if(df2 < dist_l ){  
    kdNode* nextN = nnearest(other, v,best1, dist_l, depth, pq); 
  if (nextN != NULL){
    distsq =distance(nextN->p, v); 
    if (distsq < dist_l){ 
      dist_l = distsq;
      best1 = nextN; 
      if (dist_l < pq.maxdist) 
        pq.push_node(std::make_pair(best1, dist_l)); 
    }
  }
}
return best1; 
}
