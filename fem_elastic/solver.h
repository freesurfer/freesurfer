
#ifndef H_SOLVER_H
#define H_SOLVER_H

#include <assert.h>

#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include "petscksp.h"

#include "mesh.h"
#include "cstats.h"

//----------------------------------------------------
//
// class declaration
//
//----------------------------------------------------

template<class Cstr, int n>
struct BC
{
  typedef TCoords<double,n> tCoords;
  typedef TMesh<Cstr,n> tMesh;

  tCoords pt;
  tCoords delta;
  bool    isActive;
  BC() : pt(), delta(), isActive(false)
  {}
  BC(const tCoords& _pt, const tCoords& _delta) : pt(_pt), delta(_delta),
      isActive(false)
  {}
  virtual ~BC()
  {};
  virtual bool find_candidate(tMesh* pmesh)
  {
    return false;
  }
};

template<class Cstr, int n>
struct BCNatural : public BC<Cstr,n>
{
  typedef typename BC<Cstr, n>::tCoords tCoords;
  typedef typename BC<Cstr,n>::tMesh tMesh;
  typedef TNode<n> tNode;

  tNode* pnode;
  BCNatural() : BC<Cstr, n>(), pnode(NULL)
  {}
  BCNatural(const tCoords& _pt, const tCoords& _delta) : BC<Cstr,n>(_pt, _delta), pnode(NULL)
  {}
  virtual ~BCNatural()
  {};
  virtual bool find_candidate(tMesh* pmesh)
  {
    pnode = pmesh->closest_node(this->pt);
    return (pnode!=NULL);
  }
};

template<class Cstr, int n>
struct BCMfc : public BC<Cstr, n>
{
  typedef typename BC<Cstr,n>::tCoords tCoords;
  typedef TElement<n> tElement;
  typedef typename BC<Cstr,n>::tMesh tMesh;

  tElement* pelt;
  BCMfc(const tCoords& _pt, const tCoords& _delta) : BC<Cstr,n>(_pt, _delta), pelt(NULL)
  {}
  virtual ~BCMfc()
  {};
  virtual bool find_candidate(tMesh* pmesh)
  {
    pelt = pmesh->element_at_point(this->pt);
    return (pelt!=NULL);
  }
};

//--------------------------------------------------------

template<class Cstr,int n>
class TSolver
{
public:
  typedef TCoords<int,n> tIntCoords;
  typedef TCoords<double,n> tCoords;
  typedef TMesh<Cstr,n> tMesh;
  typedef TNode<n> tNode;
  typedef TElement<n> tElement;

  typedef BC<Cstr,n> tBC;
  typedef BCNatural<Cstr,n> tBCNatural;
  typedef BCMfc<Cstr,n> tBCMfc;

  TSolver();
  TSolver(tIntCoords&);
  virtual ~TSolver();

  void clone(const TSolver&);
  void add_bc_natural(const tCoords& pt,
                      const tCoords& delta);
  void add_bc_natural(tNode* pnode,
                      const  tCoords& delta);
  void add_bc_mfc(const tCoords& pt,
                  const tCoords& delta);

  virtual int solve();
  void set_mesh(tMesh* pmesh)
  {
    m_pmesh = pmesh;
  }
  const tMesh* get_mesh() const
  {
    return m_pmesh;
  }

  void set_displayLevel(int level)
  {
    m_displayLevel = level;
  }

  int get_displayLevel() const
  {
    return m_displayLevel;
  }

  void setThreshold(double dval)
  {
    m_bcThreshold = dval;
    m_useThreshold=true;
  }

  typedef std::vector<tBC*> BcContainerType;
  typedef typename BcContainerType::const_iterator BcContainerConstIterator;
  unsigned int getBcIterators(BcContainerConstIterator& begin,
                              BcContainerConstIterator& end) const
  {
    begin = m_vBc.begin();
    end = m_vBc.end();
    return m_vBc.size();
  }

  double m_mfcWeight;

  int check_bc_error(double& dRemainingRatio);

  // mainly for debugging purposes
  // when assigning BC MFC - keep the information about
  // the element available for later probing
  typedef typename
  std::map<unsigned int, std::pair<tCoords, double> > BcMfcInfoType;
  BcMfcInfoType m_mfcInfo;


protected:
  BcContainerType m_vBc;

  Mat   m_stiffness;
  Vec   m_load;
  Vec   m_delta;

  tMesh* m_pmesh;

  int  done_bc_natural(); // distribute the BC
  int  done_bc_mfc();
  int  setup_matrix(bool showInfo=false); // assembly the stiffness matrix
  int  setup_load(); // assembly force load by reduction of the LHS
  int  setup_load_sym(); // assembly force load by
  // conditioning both rows and cols
  int  setup_load_mfc();
  int  comm_solution(); // set sol values in respective nodes

  int  add_elt_matrix(const tElement* pelt);

  int  add_elt_mfc_lhs(const tElement* pelt, tCoords& pt);
  int  add_elt_mfc_rhs(const tElement* pelt, tCoords& pt, tCoords& delta);

  int  m_displayLevel; // 0=critical, 1=important, 2=detailed

  bool m_useThreshold; // sets whether a threshold should
  // be used or not when setting the BC
  double m_bcThreshold;
};


//--------------------------------
//
// Derived class
// uses a direct solver instead of KSP
//
//  only handles natural BCs
//---------------------------------

template<class Cstr, int n>
class TDirectSolver : public TSolver<Cstr, n>
{
public:
  TDirectSolver();
  ~TDirectSolver();

  int solve();
};

//--------------------------------------------------------------------
//
// class implementation
//
//--------------------------------------------------------------------

template<class Cstr,int n>
TSolver<Cstr,n>::TSolver()
{
  m_pmesh = NULL;

  m_stiffness = 0;
  m_load = 0;
  m_delta = 0;

  m_displayLevel = 1;

  m_useThreshold = false;
  m_bcThreshold = 0.0;

  m_mfcWeight = 1.0;
}

template<class Cstr, int n>
TSolver<Cstr,n>::~TSolver()
{
  for ( typename BcContainerType::iterator it = m_vBc.begin();
        it != m_vBc.end(); ++it )
    delete *it;
  m_vBc.clear();
}

template<class Cstr,int n>
void
TSolver<Cstr,n>::clone(const TSolver& s)
{
  m_vBc   = s.m_vBc;
  //m_ticks = s.m_ticks;
  m_pmesh  = s.m_pmesh;

  m_stiffness = 0;
  m_load = 0;
  m_delta = 0;
}

template<class Cstr,int n>
void
TSolver<Cstr,n>::add_bc_natural(const tCoords& pt,
                                const tCoords& delta)
{
  tBCNatural *bc = new tBCNatural(pt,delta);
  m_vBc.push_back(bc);
}

template<class Cstr, int n>
void
TSolver<Cstr,n>::add_bc_natural(tNode* pnode,
                                const tCoords& delta)
{
  tBCNatural *bc = new tBCNatural(pnode->coords(),
                                  delta );
  bc->pnode = pnode;
  bc->pt = pnode->coords();
  m_vBc.push_back(bc);
}

template<class Cstr, int n>
void
TSolver<Cstr,n>::add_bc_mfc(const tCoords& pt,
                            const tCoords& delta)
{
  tBCMfc *bc = new tBCMfc(pt, delta);
  m_vBc.push_back(bc);
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::solve"
template<class Cstr,int n>
int
TSolver<Cstr,n>::solve()
{
  PetscErrorCode ierr;
  PetscTruth petscFlag;
  bool femPrint = false;

  ierr = PetscOptionsGetReal( NULL, "-penalty_weight",
                              &m_mfcWeight, NULL);
  CHKERRQ(ierr);
  std::cout << " penalty_weight = " << m_mfcWeight << std::endl;

  done_bc_natural();
  done_bc_mfc();

  setup_matrix(true); // display information by default for the major system

  // intermediate assembly point
  ierr = MatAssemblyBegin(m_stiffness, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_stiffness, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  {
    const unsigned int maxLen = 256;
    char buffer[maxLen];
    ierr = PetscOptionsGetString(NULL, "-fem_print",
                                 buffer, maxLen, &petscFlag);
    if ( petscFlag )
    {
      femPrint = true;
      PetscViewer viewer;
      PetscViewerBinaryOpen( PETSC_COMM_SELF, buffer,
                             FILE_MODE_WRITE, &viewer);
      MatView(m_stiffness, viewer);
      PetscViewerDestroy(viewer);
    }
  }

  setup_load_sym(); // pin down Natural conditions
  setup_load_mfc(); // setup MFC conditions

  if ( femPrint )
  {
    PetscViewer viewer;
    PetscViewerBinaryOpen( PETSC_COMM_SELF, "final_stif.bin",
                           FILE_MODE_WRITE, &viewer);
    MatView(m_stiffness, viewer);
    PetscViewerDestroy(viewer);
  }

  // final assembly point
  ierr = MatAssemblyBegin(m_stiffness, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_stiffness, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  ierr = VecAssemblyBegin(m_load);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(m_load);
  CHKERRQ(ierr);

  ierr = PetscOptionsHasName( NULL,  "-init_guess_nonzero", &petscFlag );
  CHKERRQ(ierr);
  bool initGuessNonZero = static_cast<bool>( petscFlag );
  if (initGuessNonZero )
    std::cout << " initial guess nonzero for system solving\n";

  // solve the linear system
  ierr = VecDuplicate(m_load, &m_delta);
  CHKERRQ(ierr);
  if ( initGuessNonZero )
    ierr = VecCopy(m_load, m_delta);
  CHKERRQ(ierr);

  //unused: PC pc;
  KSP ksp;
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, m_stiffness, m_stiffness, SAME_NONZERO_PATTERN);
  CHKERRQ(ierr);
  //ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  //ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp, 1.0e-9,
                          PETSC_DEFAULT, PETSC_DEFAULT,
                          PETSC_DEFAULT);
  CHKERRQ(ierr);
  if (initGuessNonZero )
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  CHKERRQ(ierr);

  // set runtime options
  ierr = KSPSetFromOptions(ksp);
  CHKERRQ(ierr);
  // solve
  ierr = KSPSolve(ksp, m_load, m_delta);
  CHKERRQ(ierr);

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp, &reason);
  CHKERRQ(ierr);
  std::cout << " convergence reason = " << reason << std::endl;

  // view solver info
  ierr = KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ(ierr);

  // check the error
  PetscInt its;
  ierr = KSPGetIterationNumber(ksp,&its);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, " Iterations %D\n", its);
  CHKERRQ(ierr);

  Vec vcheck;
  ierr = VecDuplicate(m_load, &vcheck);
  CHKERRQ(ierr);
  ierr = MatMult( m_stiffness, m_delta, vcheck);
  CHKERRQ(ierr);
  PetscScalar neg_one = -1.0;
  PetscReal norm;
  ierr = VecAXPY( vcheck, neg_one, m_load);
  CHKERRQ(ierr);
  ierr = VecNorm(vcheck, NORM_2, &norm);
  CHKERRQ(ierr);
  ierr = VecDestroy(vcheck);
  CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Absolute-Norm of error = %A\n", norm);
  CHKERRQ(ierr);

  //--------
  comm_solution(); // distribute obtained displacements to nodes

  double dremainingRatio;
  check_bc_error(dremainingRatio);
  if ( dremainingRatio > .1 && 0)
  {
    PetscViewer viewer;
    static int outCount = 0;
    char* pchBuf = new char[100];

    sprintf(pchBuf, "stif_%d.bin", outCount);
    PetscViewerBinaryOpen( PETSC_COMM_SELF,
                           pchBuf,
                           FILE_MODE_WRITE,
                           &viewer );
    MatView(m_stiffness, viewer);
    PetscViewerDestroy(viewer);

    sprintf(pchBuf, "load_%d.bin", outCount);
    PetscViewerBinaryOpen( PETSC_COMM_SELF,
                           pchBuf,
                           FILE_MODE_WRITE,
                           &viewer);
    VecView(m_load, viewer);
    PetscViewerDestroy(viewer);

    sprintf(pchBuf, "sol_%d.bin", outCount);
    PetscViewerBinaryOpen( PETSC_COMM_SELF,
                           pchBuf,
                           FILE_MODE_WRITE,
                           &viewer);
    VecView(m_delta, viewer);
    PetscViewerDestroy(viewer);

    delete[] pchBuf;
  }

  // release Petsc resources
  ierr = VecDestroy(m_delta);
  CHKERRQ(ierr);
  m_delta=0;
  ierr = VecDestroy(m_load);
  CHKERRQ(ierr);
  m_load=0;
  ierr = MatDestroy(m_stiffness);
  CHKERRQ(ierr);
  m_stiffness=0;
  ierr = KSPDestroy(ksp);
  CHKERRQ(ierr);
  ksp = 0;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::done_bc_natural"
template<class Cstr,int n>
int
TSolver<Cstr,n>::done_bc_natural()
{
  std::vector<tCoords> vdelta; // holds the point-wise diff for each BC
  bool bFailed = false;

  tNode* pnode = NULL;
  for ( typename BcContainerType::iterator it = m_vBc.begin();
        it != m_vBc.end();
        ++it )
  {
    if ( tBCNatural* bc = dynamic_cast<tBCNatural*>( *it ) )
    {
      // if node was not previously specified, find closest now
      if (!bc->pnode)
      {
        pnode = m_pmesh->closest_node(bc->pt);
        bc->pnode = pnode;
      }
      else
        pnode = bc->pnode;

      if ( pnode )
      {
        if ( !m_useThreshold ||
             (pnode->coords()- bc->pt).norm() < m_bcThreshold )
        {
          pnode->set_bc(bc->delta);
          vdelta.push_back( pnode->coords() - bc->pt );

          bc->isActive = true;

          if ( m_displayLevel>1 )
            std::cout << "setting bc " << pnode->coords()
            << " -> " << bc->delta << "\n"
            << "\t instead " << bc->pt << " -> norm = " << vdelta.back().norm()
            << std::endl;
        }
      }
      else
      {
        std::cerr
        << "TSolver::done_bc_natural -> failed to find node close to "
        << bc->pt << std::endl;
        bFailed = true;
      }
    }
  }

  if ( m_displayLevel &&
       bFailed ) std::cout << " !!!!! There were FAILED BCs\n";
  if ( m_displayLevel )
  {
    std::cout
    <<  " computing statistics for the displacement application error\n";
    double dAvgNorm = 0.0;
    for ( typename std::vector< tCoords>::const_iterator cit = vdelta.begin();
          cit != vdelta.end();
          ++cit )
      dAvgNorm += cit->norm();
    dAvgNorm /= (double)vdelta.size();
    std::cout
    << " average norm of error in placement = " << dAvgNorm << std::endl;
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::done_bc_mfc"
template<class Cstr, int n>
int
TSolver<Cstr,n>::done_bc_mfc()
{

  bool bFailed = false;

  typedef std::map< int, std::vector<int> > MfcCandidateType;
  MfcCandidateType candidates;
  MfcCandidateType::iterator mapIter;

  tElement* pelt = NULL;
  int index = 0;
  std::cout << " iterating\n";
  for ( typename BcContainerType::iterator it = m_vBc.begin();
        it != m_vBc.end(); ++it, ++index )
  {
    if ( tBCMfc* bc = dynamic_cast<tBCMfc*>(*it) )
    {
      pelt = m_pmesh->element_at_point(bc->pt);
      if ( pelt )
      {
        mapIter = candidates.find( pelt->get_id() );
        if ( mapIter == candidates.end() )
        {
          std::vector<int> vbuf;
          vbuf.push_back( index );
          candidates[ pelt->get_id() ] = vbuf;
        }
        else
          mapIter->second.push_back( index );
      }
      else
      {
        std::cerr << " TSolver::done_bc_mfc -> failed to find elt for coords "
        << bc->pt << std::endl;
        bFailed = true;
      }
    }
  } // next it
  std::cout << " done with candidates\n";

  // go through the assignment map and compute the 3D variances
  // per element
  tCoords mean;
  int active = 0;
  double dvarcova;
  for ( mapIter = candidates.begin();
        mapIter != candidates.end();
        ++mapIter )
  {
    std::vector<tCoords> vdelta;
    for ( std::vector<int>::const_iterator cit = mapIter->second.begin();
          cit != mapIter->second.end();
          ++cit )
    {
      vdelta.push_back( m_vBc[*cit]->delta );
    } // next cit
    // compute mean and variance per elt
    // return the norm of the covariance-matrix
    dvarcova = coords_statistics( vdelta, mean);
    m_mfcInfo[ mapIter->first ] = std::make_pair( mean, dvarcova );

    // get closest BC to the mean
    std::vector<int>::const_iterator citArgmin = mapIter->second.begin();
    double dMinDist(1000);
    double dCrtDist;

#if 0
    if ( dvarcova > 1.0 ) continue;
#endif

    for ( std::vector<int>::const_iterator cit = mapIter->second.begin();
          cit != mapIter->second.end();
          ++cit)
    {
      // need to write a routine to invert the covariance
      // matrix 3x3 - should be direct
      // use determinants, i guess
      dCrtDist = ( m_vBc[*cit]->delta - mean).norm();
      if ( dCrtDist < dMinDist )
      {
        dMinDist = dCrtDist;
        citArgmin = cit;
      }
    } // next cit

    // assign BC
    ++active;
    m_vBc[*citArgmin]->isActive = true;
    dynamic_cast<tBCMfc*>(m_vBc[*citArgmin])->pelt =
      m_pmesh->fetch_elt(mapIter->first);
  } // next mapIter

  std::cout << " Active BCs = " << active << std::endl
  << " Total BCs = " << m_vBc.size() << std::endl;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::setup_matrix"
template<class Cstr,int n>
int
TSolver<Cstr,n>::setup_matrix(bool showInfo)
{
  PetscErrorCode ierr;
  int n_eqs = m_pmesh->get_no_nodes() * n; // template par = dim
  if (showInfo)  std::cout << " no-eqs = " << n_eqs << std::endl;
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n_eqs,n_eqs,
                         100, NULL,
                         &m_stiffness);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(m_stiffness);
  CHKERRQ(ierr);
  ierr = MatSetOption(m_stiffness, MAT_SYMMETRIC);
  CHKERRQ(ierr);
  ierr = MatSetOption(m_stiffness, MAT_SYMMETRY_ETERNAL);
  CHKERRQ(ierr);

  int iold_val = -1;
  for ( size_t i= size_t(0);
        i<m_pmesh->get_no_elts();
        ++i )
  {
    if ( showInfo )
    {
      int idone = (int)floor( double(i) / m_pmesh->get_no_elts() * 100.0 );
      if ( idone%5==0 && idone!=iold_val )
      {
        iold_val=idone;
        std::cout << "\t percentage done= " << idone << std::endl;
      }
    }
    add_elt_matrix(m_pmesh->get_elt(i));
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::add_elt_matrix"
template<class Cstr, int n>
int
TSolver<Cstr,n>::add_elt_matrix(const tElement* pelt)
{
  PetscErrorCode ierr;
  int *id = new int[ pelt->no_nodes() ];
  tNode* pnode = NULL;

  for (int i= 0;
       i<pelt->no_nodes();
       ++i )
  {
    if (!pelt->get_node(i,&pnode))
    {
      std::cerr << "TSolver::add_elt_matrix -> err 1\n";
      exit(1);
    }
    id[i] = pnode->get_id();
  }

  SmallMatrix elt_matrix = pelt->get_matrix();

  PetscScalar* values = new PetscScalar[ elt_matrix.rows()
                                         * elt_matrix.cols() ];
  int *pindex_1 = new int[ elt_matrix.rows() ];
  int *pindex_2 = new int[ elt_matrix.cols() ];

  int a,b;
  for (int i= 0; i<pelt->no_nodes(); ++i)
    for (int j= 0; j<pelt->no_nodes(); ++j)
      for (int k=0; k<n; ++k)
        for (int l=0; l<n; ++l)
        {
          a = n*i + k;
          b = n*j + l;

          assert(id[i]>=0);
          assert(id[j]>=0);

          pindex_1[a] = n*id[i]+k;
          pindex_2[b] = n*id[j]+l;

          values[ a + n* pelt->no_nodes() *b ] = (float)elt_matrix(a,b);
        }

  ierr = MatSetValues( m_stiffness,
                       elt_matrix.rows(), pindex_1,
                       elt_matrix.cols(), pindex_2,
                       values,
                       ADD_VALUES);
  CHKERRQ(ierr);

  delete[] id;
  delete[] values;
  delete[] pindex_1;
  delete[] pindex_2;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::add_elt_mfc_lhs"
template<class Cstr, int n>
int
TSolver<Cstr,n>::add_elt_mfc_lhs(const tElement* pelt,
                                 tCoords& pt)
{
  PetscErrorCode ierr;
  int *id = new int[ pelt->no_nodes() ];
  tNode* pnode = NULL;

  for (int i=0; i<pelt->no_nodes(); ++i)
  {
    if ( !pelt->get_node(i,&pnode) )
    {
      std::cerr << " TSolver::add_elt_mfc -> err 2\n";
      exit(1);
    }
    id[i] = pnode->get_id();
  }

  SmallMatrix elt_matrix(n, n*pelt->no_nodes());

  SmallMatrix bufMatrix;
  for ( int i=0; i<pelt->no_nodes(); ++i)
  {
    bufMatrix.identity(n);
    bufMatrix *= pelt->shape_fct(i, pt);

    elt_matrix.set_block( bufMatrix,
                          0, i*n );
  }

  // obtain the matrix by LEFT multiplying with the transpose
  bufMatrix = elt_matrix.transpose() * elt_matrix;
  bufMatrix *= m_mfcWeight;

  int *pindex_1 = new int[ bufMatrix.rows() ];
  int *pindex_2 = new int[ bufMatrix.cols() ];

  PetscScalar* values = new PetscScalar[ bufMatrix.rows() *
                                         bufMatrix.cols() ];

  int a,b;
  for (int i=0; i<pelt->no_nodes(); ++i)
    for (int j=0; j<pelt->no_nodes(); ++j)
      for (int k=0; k<n; ++k)
        for (int l=0; l<n; ++l)
        {
          a = n*i + k;
          b = n*j + l;

          pindex_1[a] = n*id[i] + k;
          pindex_2[b] = n*id[j] + l;

          values[ a + n* pelt->no_nodes() * b] = bufMatrix(a,b);
        } // next i,j,k,l


  ierr = MatSetValues( m_stiffness,
                       bufMatrix.rows(), pindex_1,
                       bufMatrix.cols(), pindex_2,
                       values,
                       ADD_VALUES);
  CHKERRQ(ierr);

  delete[] id;
  delete[] pindex_1;
  delete[] pindex_2;
  delete[] values;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::setup_load_mfc"
template<class Cstr, int n>
int
TSolver<Cstr,n>::setup_load_mfc()
{
  for ( typename BcContainerType::iterator it = m_vBc.begin();
        it != m_vBc.end(); ++it )
  {
    if ( !(*it)->isActive ) continue;

    if ( tBCMfc* bc = dynamic_cast<tBCMfc*>( &(*(*it)) ) )
    {
      add_elt_mfc_lhs( bc->pelt, bc->pt );
      add_elt_mfc_rhs( bc->pelt, bc->pt, bc->delta);
    }
  } // next it

  std::cout << " after setup_load_mfc\n";
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::add_elt_mfc_rhs"
template<class Cstr, int n>
int
TSolver<Cstr, n>::add_elt_mfc_rhs(const tElement* pelt,
                                  tCoords& pt,
                                  tCoords& delta)
{
  PetscErrorCode ierr;
  tNode* pnode = NULL;

  SmallMatrix elt_matrix(n, n* pelt->no_nodes() );
  SmallMatrix bufMatrix;

  for (int i=0; i<pelt->no_nodes(); ++i)
  {
    bufMatrix.identity(n);
    bufMatrix *= pelt->shape_fct(i, pt);

    elt_matrix.set_block( bufMatrix,
                          0, i*n);
  } // next i

  elt_matrix.inplace_transpose();

  SmallMatrix mfc_rhs(n, 1);
  for (int i=0; i<n; ++i)
    mfc_rhs(i,0) = delta(i) * m_mfcWeight;

  bufMatrix = elt_matrix * mfc_rhs;

  std::map<int, double> mrhs;

  for (int i=0; i<pelt->no_nodes(); ++i)
  {
    if ( !pelt->get_node(i, &pnode) )
    {
      std::cerr << " TSolver::add_elt_mfc_rhs -> err 3\n";
      exit(1);
    }
    for (int j=0; j<n; ++j)
      mrhs[ n* pnode->get_id() + j ] = bufMatrix( n*i + j, 0);
  } // next i

  int    *indices = new int[ mrhs.size() ];
  double *values  = new double[ mrhs.size() ];

  int index =0;
  for ( typename std::map<int, double>::const_iterator cit = mrhs.begin();
        cit != mrhs.end(); ++cit, ++index )
  {
    indices[ index ] = cit->first;
    values[ index ] = cit->second;
  } // next cit

  ierr = VecSetValues( m_load, (int)mrhs.size(),
                       indices, values,
                       ADD_VALUES);
  CHKERRQ(ierr);

  delete[] indices;
  delete[] values;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::comm_solution"
template<class Cstr,int n>
int
TSolver<Cstr,n>::comm_solution()
{
  PetscErrorCode ierr;

  int no_eqs = m_pmesh->get_no_nodes() *n;

  PetscScalar *pdelta = new PetscScalar[no_eqs];
  ierr = VecGetArray(m_delta, &pdelta);
  CHKERRQ(ierr);

  tNode* pnode = NULL;
  int count = 0;
  for (size_t i=size_t(0); i<m_pmesh->get_no_nodes(); ++i,count++)
  {
    pnode = NULL;
    m_pmesh->get_node(i,&pnode);
    if ( !pnode )
    {
      std::cerr << "TSolver::commm_solution -> err\n";
      exit(1);
    }
    for (int j=0; j<n; ++j)
      pnode->set_dof_val(j, pdelta[ pnode->get_id()*n + j]);
  }
  ierr = VecRestoreArray(m_delta, &pdelta);
  CHKERRQ(ierr);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::setup_load"
template<class Cstr,int n>
int
TSolver<Cstr,n>::setup_load()
{
  PetscErrorCode ierr;
  tNode* pnode;
  int no_eqs = n * m_pmesh->get_no_nodes();

  ierr = MatSetOption(m_stiffness, MAT_NO_NEW_NONZERO_LOCATIONS);
  CHKERRQ(ierr);

  std::map<int, double> mrhs;

  for (size_t i=size_t(0); i<m_pmesh->get_no_nodes(); ++i)
  {
    if ( !m_pmesh->get_node(i, &pnode) )
    {
      std::cerr << "TSolver::setup_load -> requested node out of range\n";
      exit(1);
    }

    for (int j=0; j<n; ++j)
    {
      if (pnode->is_dof_active(j))
        mrhs[ n* pnode->get_id() + j ] = pnode->get_dof(j);
    }
  }

  // create the RHS vector
  // by default, all the values will be NULL
  ierr = VecCreate( PETSC_COMM_SELF, &m_load);
  CHKERRQ(ierr);
  ierr = VecSetSizes(m_load, PETSC_DECIDE, no_eqs);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(m_load);
  CHKERRQ(ierr);

  int* indices = new int[ mrhs.size()];
  double* values = new double[ mrhs.size()];

  int i=0;
  for ( typename std::map<int,double>::const_iterator
        cit = mrhs.begin();
        cit!= mrhs.end();
        ++cit ,++i)
  {
    indices[i] = cit->first;
    values[i] = cit->second;
  }

  // create index set

  IS is;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, (int)mrhs.size(),
                         indices, &is);
  CHKERRQ(ierr);
  ierr = MatZeroRowsIS(m_stiffness, is, 1.0);
  CHKERRQ(ierr);
  ierr = ISDestroy(is);
  CHKERRQ(ierr);

  ierr = VecSetValues(m_load, (int)mrhs.size(),
                      indices, values,
                      INSERT_VALUES);
  CHKERRQ(ierr);

  ierr = VecAssemblyBegin(m_load);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(m_load);
  CHKERRQ(ierr);

  delete[] indices;
  delete[] values;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::setup_load_sym"
template <class Cstr, int n>
int
TSolver<Cstr,n>::setup_load_sym()
{
  PetscErrorCode ierr;
  //unused: tNode* pnode;

  int no_eqs = n * m_pmesh->get_no_nodes();

  ierr = MatSetOption(m_stiffness, MAT_NO_NEW_NONZERO_LOCATIONS);
  CHKERRQ(ierr);
  std::map<int,double> mrhs;

#if 0
  for (size_t i=size_t(0); i<m_pmesh->get_no_nodes(); ++i)
  {
    if ( !m_pmesh->get_node(i, &pnode) )
    {
      std::cerr << "TSolver::setup_load_sym -> requested node out of range\n";
      exit(1);
    }

    for (int j=0; j<n; ++j)
    {
      if (pnode->is_dof_active(j))
        mrhs[ n* pnode->get_id() + j ] = pnode->get_dof(j);
    }      // next j
  } // next i
#endif
  if ( m_displayLevel>1) std::cout << " done_bc_natural size of container = "
    << m_vBc.size() << std::endl;
  ;
  for (typename BcContainerType::iterator it = m_vBc.begin();
       it != m_vBc.end(); ++it )
  {
    if ( !(*it)->isActive ) continue;

    if ( tBCNatural* bc = dynamic_cast<tBCNatural*>(*it) )
    {
      for (int j=0; j<n; ++j)
        mrhs[ n* bc->pnode->get_id() + j ] = bc->pnode->get_dof(j);
    }
  }
  if ( m_displayLevel>1) std::cout << " after building the map\n";

  std::cout << " LOAD size = " << mrhs.size() << std::endl;

  // create RHS vector
  ierr = VecCreate( PETSC_COMM_SELF, &m_load);
  CHKERRQ(ierr);
  ierr = VecSetSizes(m_load, PETSC_DECIDE, no_eqs);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(m_load);
  CHKERRQ(ierr);

  int *indices = new int[ mrhs.size() ];
  double *values = new double[ mrhs.size() ];

  int i=0;
  for ( typename std::map<int,double>::const_iterator
        cit = mrhs.begin();
        cit!= mrhs.end();
        ++cit ,++i)
  {
    indices[i] = cit->first;
    values[i] = cit->second;
  }

  // compute RHS vector
  Vec vecBcs;
  ierr = VecCreate( PETSC_COMM_SELF, &vecBcs);
  CHKERRQ(ierr);
  ierr = VecSetSizes( vecBcs, PETSC_DECIDE, no_eqs);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(vecBcs);
  CHKERRQ(ierr);

  PetscScalar* rhsValues = new PetscScalar[ mrhs.size()];

  for ( i=0; i<(int)mrhs.size(); ++i)
    rhsValues[i] = -values[i];

  ierr = VecSetValues( vecBcs, (int)mrhs.size(), indices,
                       rhsValues, INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(vecBcs);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vecBcs);
  CHKERRQ(ierr);

  ierr = MatMult( m_stiffness, vecBcs, m_load);
  CHKERRQ(ierr);
  ierr = VecDestroy(vecBcs);
  CHKERRQ(ierr);
  ierr = VecSetValues( m_load, (int)mrhs.size(),
                       indices, values, INSERT_VALUES);
  CHKERRQ(ierr);

  // condition the matrix
  IS is;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, (int)mrhs.size(),
                         indices, &is);
  CHKERRQ(ierr);
  ierr = MatZeroRowsIS(m_stiffness, is, 1.0);
  CHKERRQ(ierr);
  ierr = MatTranspose( m_stiffness, PETSC_NULL);
  CHKERRQ(ierr);
  ierr = MatZeroRowsIS(m_stiffness, is, 1.0);
  CHKERRQ(ierr);

  ierr = ISDestroy(is);
  CHKERRQ(ierr);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::check_bc_error"
template <class Cstr, int n>
int
TSolver<Cstr,n>::check_bc_error(double& dRemainingRatio)
{
  // go through the BCs and measure the error
  double dSum(.0), dSumConditional(.0), dSumInitial(.0), dval;
  int count(0), countConditional(0);

  tCoords img;
  unsigned int countInvalid = 0;
  for ( typename BcContainerType::const_iterator cit = m_vBc.begin();
        cit != m_vBc.end(); ++cit )
  {
    img = m_pmesh->dir_img( (*cit)->pt );
    if ( !img.isValid() )
    {
      ++countInvalid;
      continue;
    }
    // if a topology problem is observed, no point carrying on

    dval = ( img - (*cit)->pt - (*cit)->delta ).norm();

    dSumInitial += (*cit)->delta.norm();

    dSum += dval;
    ++count;
    if ( (*cit)->isActive )
    {
      dSumConditional += dval;
      countConditional++;
    }
  } // next cit
  std::cout << " countInvalid = " << countInvalid
  << " general-count = " << count << std::endl;

  if ( count )
  {
    std::cout << " Average of the error norm = "
    << dSum /(double)count << std::endl
    << " Initial error = " << dSumInitial / (double)count << std::endl;
  }
  else
    std::cout << " count = 0 !?!\n";

  if ( countConditional )
    std::cout << " Average of the error norm conditional = "
    << dSumConditional / (double)countConditional << std::endl;
  else
    std::cout << " countConditional = 0 !?!?\n";

  dRemainingRatio = dSum / dSumInitial;

  return 0;

}

//-----------------------------------------------------------------------
//
//

template<class Cstr, int n>
TDirectSolver<Cstr,n>::TDirectSolver()
    : TSolver<Cstr,n>()
{
  this->m_displayLevel = 1;
}

template<class Cstr, int n>
TDirectSolver<Cstr,n>::~TDirectSolver()
{}

#undef __FUNCT__
#define __FUNCT__ "TDirectSolver::solve"
template<class Cstr, int n>
int
TDirectSolver<Cstr,n>::solve()
{
  PetscErrorCode ierr;

  this->done_bc_natural();

  this->setup_matrix();

  // intermediate assembly point
  ierr = MatAssemblyBegin(this->m_stiffness, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(this->m_stiffness, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  this->setup_load_sym();

  // final assembly point
  ierr = MatAssemblyBegin(this->m_stiffness, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(this->m_stiffness, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  ierr = VecAssemblyBegin(this->m_load);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(this->m_load);
  CHKERRQ(ierr);

  // solve linear system
  ierr = VecDuplicate(this->m_load, &this->m_delta);
  CHKERRQ(ierr);
  ierr = VecCopy(this->m_load, this->m_delta);
  CHKERRQ(ierr);

  PC pc;
  KSP ksp;
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,
                         this->m_stiffness,
                         this->m_stiffness,
                         SAME_NONZERO_PATTERN);
  CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc);
  CHKERRQ(ierr);
  ierr = PCSetType(pc, PCJACOBI);
  CHKERRQ(ierr);
  ierr = PCFactorSetUseInPlace(pc);
  ierr = KSPSetType(ksp, KSPCG);
  CHKERRQ(ierr);

  ierr = KSPSetOptionsPrefix(ksp, "direct_");
  ierr = KSPSetFromOptions(ksp);
  CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = KSPSolve(ksp, this->m_load, this->m_delta);
  CHKERRQ(ierr);

  {
    static std::map<int,std::string> kspConvergenceReason;
    kspConvergenceReason[KSP_CONVERGED_RTOL] = "ksp-converged-rtol";
    kspConvergenceReason[KSP_CONVERGED_ATOL] = "ksp-converged-atol";
    kspConvergenceReason[KSP_CONVERGED_ITS]  = "ksp-converged-its";
    kspConvergenceReason[KSP_CONVERGED_CG_NEG_CURVE] =
      "ksp-converged-stcg-neg-curve";
    kspConvergenceReason[KSP_CONVERGED_CG_CONSTRAINED] =
      "ksp-converged-stcg-constrained";
    kspConvergenceReason[KSP_CONVERGED_STEP_LENGTH] =
      "ksp-converged-step-length";
    kspConvergenceReason[KSP_CONVERGED_HAPPY_BREAKDOWN] =
      "ksp-converged-happy-breakdown";
    kspConvergenceReason[KSP_DIVERGED_NULL] = "ksp-diverged-null";
    kspConvergenceReason[KSP_DIVERGED_ITS] = "ksp-diverged-its";
    kspConvergenceReason[KSP_DIVERGED_DTOL] = "ksp-diverged-dtol";
    kspConvergenceReason[KSP_DIVERGED_BREAKDOWN] = "ksp-diverged-breakdown";
    kspConvergenceReason[KSP_DIVERGED_BREAKDOWN_BICG]=
      "ksp-diverged-breakdown-bicg";
    kspConvergenceReason[KSP_DIVERGED_NONSYMMETRIC] =
      "ksp-diverged-nonsymmetric";
    kspConvergenceReason[KSP_DIVERGED_INDEFINITE_PC]=
      "ksp-diverged-indefinite-pc";
    kspConvergenceReason[KSP_DIVERGED_NAN]="ksp-diverged-nan";
    kspConvergenceReason[KSP_DIVERGED_INDEFINITE_MAT]=
      "ksp-diverged-indefinite-mat";
    kspConvergenceReason[KSP_CONVERGED_ITERATING]="ksp-converged-iterating";

    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(ksp, &reason);
    CHKERRQ(ierr);
    if (this->m_displayLevel)
      std::cout << "DirectSolverConvergence = "
      << kspConvergenceReason[reason] << std::endl;
  }

  // distribute obtained displacements to nodes
  this->comm_solution();

  // release Petsc resources
  ierr = VecDestroy(this->m_delta);
  CHKERRQ(ierr);
  this->m_delta = 0;
  ierr = VecDestroy(this->m_load);
  CHKERRQ(ierr);
  this->m_load=0;
  ierr = MatDestroy(this->m_stiffness);
  CHKERRQ(ierr);
  this->m_stiffness=0;
  ierr = KSPDestroy(ksp);
  CHKERRQ(ierr);
  ksp = 0;

  return 0;
}

#endif
