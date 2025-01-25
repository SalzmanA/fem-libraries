#define MAX_REFINE 0
#define MAX_DAM 1.
//#define MAX_DAM 0.
#define USE_VOLUME
// bugged in 0.8.0 in pack_coefficient_entity #define USE_SURF
//#define USE_SQUARE
#define USE_TRAC
//#define PROF_KERNEL
#define USE_BOOMERAMG
#define USE_BOOMERAMG_FINE_TUNNING 1
//#define IN_COMP

#define USE_ECST

#ifdef HAS_ADIOS2
#define USE_ADIOS_FOR_OUTPUT
#else
#define USE_XDMF_FOR_OUTPUT
//#define USE_VTK_FOR_OUTPUT
#endif

#ifdef USE_TRAC
#ifdef USE_ECST
#undef USE_ECST
#endif
#endif

#include <unordered_set>
#include <string>
#include <dolfinx/fem/petsc.h>
#ifdef SYMB_SYM
#include "asym_symb_sym.h"
#define SPACE_F functionspace_form_asym_symb_sym_F
#define FORM_F form_asym_symb_sym_F
#define FORM_J form_asym_symb_sym_J
#define EXPR0_F expression_asym_symb_sym_0
#define EXPR1_F expression_asym_symb_sym_1
#define EXPR2_F expression_asym_symb_sym_2
#define SEXT "symb_sym";
#elif defined (SYMB)
#include "asym_symb.h"
#define SPACE_F functionspace_form_asym_symb_F
#define FORM_F form_asym_symb_F
#define FORM_J form_asym_symb_J
#define EXPR0_F expression_asym_symb_0
#define EXPR1_F expression_asym_symb_1
#define EXPR2_F expression_asym_symb_2
#define SEXT "symb";
#elif defined (MANUAL)
#include "asym_manual.h"
#define SPACE_F functionspace_form_asym_manual_F
#define FORM_F form_asym_manual_F
#define FORM_J form_asym_manual_J
#define EXPR0_F expression_asym_manual_0
#define EXPR1_F expression_asym_manual_1
#define EXPR2_F expression_asym_manual_2
#define SEXT "manual";
#elif defined (UFL_POTENTIAL)
#include "asym_ufl.h"
#define SPACE_F functionspace_form_asym_ufl_F
#define FORM_F form_asym_ufl_F
#define FORM_J form_asym_ufl_J
#define EXPR0_F expression_asym_ufl_0
#define EXPR1_F expression_asym_ufl_1
#define EXPR2_F expression_asym_ufl_2
#define SEXT "ufl";
#endif
#include "dolfinx.h"
#include "dolfinx/graph/partitioners.h"
#include "dolfinx/io/XDMFFile.h"
#include <dolfinx/io/ADIOS2Writers.h>
#include <dolfinx/nls/NewtonSolver.h>
#include "petscpctypes.h"

#define SL(i, m)                                                                                                                \
   cout << "| " << fixed << setprecision(5) << setw(12) << min_measure[i] << " | " << setw(12) << max_measure[i] << " | "       \
        << setw(12) << sqrt(ecart_type_measure[i] * nb_proc - avg_measure[i] * avg_measure[i]) / nb_proc << " | " << setw(12)   \
        << 100*sqrt(ecart_type_measure[i] * nb_proc - avg_measure[i] * avg_measure[i]) / avg_measure[0] << " | " << setprecision(5) \
        << setw(12) << avg_measure[i] / nb_proc << " | " << setprecision(1) << setw(12)                                         \
        << 100 * avg_measure[i] / avg_measure[0];                                                                               \
   cout << " | " << setw(42) << m << " |" << endl;

#ifdef PROF_KERNEL
#define NB_MEASURE 18
#else
#define NB_MEASURE 16
#endif
using namespace std;

#ifdef IN_COMP
#include <algorithm>
#include <fstream>
typedef array<double,4> quarted_t;
typedef vector<quarted_t> container_t;
#define EPS 1.e-10
bool double_eq(const double &a, const double &b, double eps=EPS)
{
   if ((a-b)>eps)return false;
   if ((b-a)>eps)return false;
   return true;
}
#endif


int main(int argc, char *argv[])
{
   //== init ==================================================
   // init logging backend of dolfinx (loguru)
   // tuned by -dolfinx_loglevel <level> option
   dolfinx::init_logging(argc, argv);
   // init petsc environment
   // imply init distributed environement
   PetscInitialize(&argc, &argv, nullptr, nullptr);
   const std::size_t nb_proc = dolfinx::MPI::size(MPI_COMM_WORLD);
   const int mpi_rank = dolfinx::MPI::rank(MPI_COMM_WORLD);
   const bool do_print = (!mpi_rank);
   double measure[NB_MEASURE];
   double start_all = MPI_Wtime();

   // change scope so that all dolfinx object (or related) get destroy before MPI finalization
   {
      double start = MPI_Wtime();
      std::string thread_name = "RANK: " + std::to_string(mpi_rank);
      loguru::set_thread_name(thread_name.c_str());
      if (do_print)
         loguru::g_stderr_verbosity =
#ifdef NDEBUG
             loguru::Verbosity_WARNING;
#else
             loguru::Verbosity_INFO;
#endif

      // init output per proc
      string no = "proc_" + std::to_string(mpi_rank) + "_output.txt";
      if (mpi_rank > 50)
         auto fd = freopen("/dev/null", "w", stdout);
      else
         auto fd = freopen(no.c_str(), "w", stdout);
      measure[1] += MPI_Wtime() - start;  // 1. Initialize

      //== types =================================================
      typedef PetscScalar scalar;
      typedef dolfinx::scalar_value_type_t<scalar> scalar_dolf;

      //== triangle element =========================================
      auto tria = dolfinx::fem::CoordinateElement<double>(dolfinx::mesh::CellType::triangle, 1);
      int dim = dolfinx::mesh::cell_dim(dolfinx::mesh::CellType::triangle);
      int dim1 = dolfinx::mesh::cell_dim(dolfinx::mesh::CellType::interval);
      int dim0 = dolfinx::mesh::cell_dim(dolfinx::mesh::CellType::point);

      //== read mesh =============================================
      start = MPI_Wtime();
#ifdef USE_SQUARE
      dolfinx::io::XDMFFile file(MPI_COMM_WORLD, "data/square.xdmf", "r");
#else
      dolfinx::io::XDMFFile file(MPI_COMM_WORLD, "data/neper_dam.xdmf", "r");
#endif
      dolfinx::mesh::Mesh<double> mesh = file.read_mesh(tria, dolfinx::mesh::GhostMode::none, "neper_dam");
      dolfinx::mesh::MeshTags<int32_t> cell_tag = file.read_meshtags(mesh, "neper_dam_cells");
      mesh.topology()->create_entities(dim1);
      dolfinx::mesh::MeshTags<int32_t> edge_tag = file.read_meshtags(mesh, "neper_dam_facets");
      auto pmesh = std::make_shared<dolfinx::mesh::Mesh<double>>(mesh);
      measure[2]+=MPI_Wtime()-start; // 2.1 Read the mesh
      start = MPI_Wtime();
      // refine
      for (short i = 0; i < MAX_REFINE; ++i) 
      {
         auto nm = dolfinx::refinement::plaza::refine(*pmesh, false,dolfinx::refinement::plaza::Option::parent_cell_and_facet);
         pmesh.reset();
         pmesh=std::make_shared<dolfinx::mesh::Mesh<double>>(std::get<0>(nm));
         auto nt = pmesh->topology();
         auto nct=dolfinx::refinement::transfer_cell_meshtag(cell_tag,*nt,std::get<1>(nm));
         cell_tag=dolfinx::mesh::MeshTags<int32_t>(nt,dim,nct[0],nct[1]);
         nt->create_entities(dim1);
         auto nft=dolfinx::refinement::transfer_facet_meshtag(edge_tag,*nt,std::get<1>(nm),std::get<2>(nm));
         edge_tag = dolfinx::mesh::MeshTags<int32_t>(nt, dim1, nft[0], nft[1]);
#ifdef USE_SQUARE
         cout<<"Refinement "<<i<<" :"<<endl;
         for (auto v : cell_tag.indices()) cout<<" "<<v;
         cout<<endl;
         for (auto v : cell_tag.values()) cout<<" "<<v;
         cout<<endl;
#endif
      }
      measure[3]+=MPI_Wtime()-start; // 2.2 Refine the mesh
      start = MPI_Wtime();
      auto topo = pmesh->topology();
      topo->create_entities(dim0);
      topo->create_connectivity(dim1, dim0);
      topo->create_connectivity(dim0, dim1);
      auto edges = topo->connectivity(dim1, dim0);
      measure[2]+=MPI_Wtime()-start; // 2.1 Read the mesh

      //== mesh distribution statistic ===========================
      {
         topo->create_connectivity(dim, dim0);
         auto cells = topo->connectivity(dim, dim0);
         cout << " local mesh statistic:" << endl;
         cout << "*" << cells->num_nodes() << " elements" << endl;
#ifdef USE_SQUARE
         for (std::int32_t l = 0; l < cells->num_nodes(); ++l)
         {
            cout << "TRI " << l << " adj ";
            auto link = cells->links(l);
            for (auto n : link) cout << " " << n;
            cout << endl;
         }
         topo->create_connectivity(dim0, dim);
         auto nodes = topo->connectivity(dim0, dim);
         cout <<"*"<<nodes->num_nodes() << " nodes" << endl;
         auto geo =pmesh->geometry().x();
         MDSPAN_IMPL_STANDARD_NAMESPACE::mdspan nxyz(&geo[0],nodes->num_nodes(),3);
         cout << "Nodes coordinates " << endl;
         for (size_t l = 0; l < nodes->num_nodes(); ++l)
         {
            cout << " " << l;
            for (size_t i=0; i < 3; ++i) cout << " " << nxyz(l,i);
            cout << endl;
         }
         cout<<"*"<< edges->num_nodes()<<" edges"<<endl;
         for (std::int32_t l = 0; l < edges->num_nodes(); ++l)
         {
            auto link = edges->links(l);
            cout << l;
            for (auto nidx : link) cout << " " << nidx;
            cout << endl;
         }
         std::string dim_n[] = {"nodes", "edges", "faces"};
         for (int dim_i=0;dim_i<3;++dim_i)
         {
            auto &ime = *topo->index_map(dim_i);
            cout << dim_n[dim_i] + " : nb ghost " << ime.num_ghosts() << " local size " << ime.size_local() << " global size "
                 << ime.size_global() << endl;
            auto lg = ime.global_indices();
            cout << dim_n[dim_i] + " global index";
            for (auto g : lg) cout << " " << g;
            cout << endl;
            auto ll = ime.local_range();
            cout << dim_n[dim_i] + " range in global indexing of local indexes";
            for (auto g : ll) cout << " " << g;
            cout << endl;
         }
#endif
      }

      //== misc parameter ============================================
      const double eps = 1.0e-8;
      const double eps1 = 1.-1.0e-8;

      //== Constant associated to the forms ==========================
      // UFC  nu (Poisson ratio) is defined in  python form 
      const std::map<std::string, std::shared_ptr<const dolfinx::fem::Constant<scalar>>> constants = {
          {"nu", std::make_shared<const dolfinx::fem::Constant<scalar>>(dolfinx::fem::Constant<scalar>(0.3))}
#ifdef USE_SURF
#ifdef USE_TRAC
          ,{"t", std::make_shared<const dolfinx::fem::Constant<scalar>>(dolfinx::fem::Constant<scalar>(10000.))}
#else
          ,{"t", std::make_shared<const dolfinx::fem::Constant<scalar>>(dolfinx::fem::Constant<scalar>(-100.))}
#endif
#endif
      };

      start = MPI_Wtime();
      //== function space associated to the problem ==================
      // UFC SPACE_F comes from c files generated from python form file with ffcx
      // u,delta_u and f share the same function space
      auto V = std::make_shared<dolfinx::fem::FunctionSpace<scalar_dolf>>(
          dolfinx::fem::create_functionspace(SPACE_F, "u", pmesh));
      // DG space used to store Young modulus E
      auto ES = std::make_shared<dolfinx::fem::FunctionSpace<scalar_dolf>>(
          dolfinx::fem::create_functionspace(SPACE_F, "E", pmesh));
      // scalar space used to store d (and durty one)
      auto DS = std::make_shared<dolfinx::fem::FunctionSpace<scalar_dolf>>(
          dolfinx::fem::create_functionspace(SPACE_F, "d", pmesh));
      measure[5]+=MPI_Wtime()-start; // 3.1 Define space

#ifdef USE_SQUARE
      {
         cout << "Values size for space V "<<V->value_size()<<endl;
         cout << "Values size for space ES "<<ES->value_size()<<endl;
         cout << "Values size for space DS "<<DS->value_size()<<endl;
         cout << "Dof coordinates for space DS/V";
         auto tdc = DS->tabulate_dof_coordinates(false);
         std::int32_t p = 0;
         for (auto nxyz : tdc)
         {
            if (!(p % 3)) cout << endl << "dof " << p / 3;
            cout << " " << nxyz;
            ++p;
         }
         cout << endl;
         auto dm = DS->dofmap();
         auto m = dm->map();
         auto &im = *dm->index_map;
         cout << "Dof id per cell for space DS/V" << endl;
         for (std::int32_t l = 0; l < m.extent(0); ++l)
         {
            cout << "cell " << l<<" dofs ";
            for (std::size_t p = 0; p < m.extent(1); ++p) cout << " " << m(l, p);
            cout << endl;
         }
         cout<<"DS : nb ghost "<<im.num_ghosts()<<" local size "<<im.size_local()<<" global size "<<im.size_global()<<" bs "<<dm->index_map_bs()<<endl;
         auto lg = im.global_indices();
         cout << "DS dof global index";
         for (auto g : lg) cout << " " << g;
         cout<<endl;
         auto ll = im.local_range();
         cout << "DS dof range in global indexing of local indexes"; 
         for (auto g : ll) cout << " " << g;
         cout<<endl;
      }
#endif

      start = MPI_Wtime();
      //== create damage d field =====================================
      auto d = std::make_shared<fem::Function<scalar>>(DS);
      d->name = "d";
      d->x()->set(0.);
      {
     

         // For ghost needs to put in place a correspondence of index between nodes and dofs as we
         // gone use dofs/nodes for computation. To simplify a full correspondence is set for local
         // and ghost. Pass by global !!! Ouch !!! Don't see for now how to do that better
         // !!! and it work only because Lagrange1 scalar dof and nodes are linked appropriately
         // using ghosts may be quicker: only ask for ghosts, not global_indices so smaller vector
         auto dm = DS->dofmap();
         auto m = dm->map();
         auto &im_space = *dm->index_map;
         auto &im_mesh = *topo->index_map(dim0);
         int32_t nloc = im_space.size_local();
         int32_t nglob = im_space.num_ghosts();
         int32_t nlg = nloc+nglob;
         assert(nloc == im_mesh.size_local());
         assert(nglob == im_mesh.num_ghosts());
         std::vector<std::int64_t> global_mesh_indices = im_mesh.global_indices();
         std::vector<std::int32_t> corresp(nlg);
         size_t p = 0;
         for (auto &idx : corresp) idx=p++;
         std::span<std::int32_t> loc(corresp.data()+nloc,nglob);
         std::span<std::int64_t> glob(global_mesh_indices.data()+nloc,nglob);
         im_space.global_to_local(glob,loc);
         global_mesh_indices.clear();

         // extra distributed vector to operate directly on mesh index
         dolfinx::la::Vector<scalar> nd_valv(topo->index_map(dim0), 1);
         auto nd_val=nd_valv.mutable_array();

         // number of edge per owned node
         std::vector<scalar> num_edges_per_nodes(nloc);

         // collect edges chosen arbitrarily
         std::vector<std::int32_t> edges_damaged;
#ifdef USE_SQUARE
         std::vector<std::int32_t> tag_edges_damaged = {4};
         //std::vector<std::int32_t> tag_edges_damaged = {3, 4, 7};
#else
         std::vector<std::int32_t> tag_edges_damaged = {148, 342, 333, 19,  380, 408, 328, 329, 325, 323,
                                                        96,  97,  531, 4,   471, 234, 235, 184, 236, 419,
                                                        350, 332, 364, 176, 77,  333, 341, 343, 144, 143};
#endif
         for (auto tag : tag_edges_damaged)
         {
            auto s = edge_tag.find(tag);
            edges_damaged.insert(edges_damaged.end(), s.begin(), s.end());
         }
         // collect nodes with damage
         std::unordered_set<std::int32_t> node_damaged;
         for (auto l : edges_damaged)
         {
            auto link = edges->links(l);
            for (auto nidx : link) 
            node_damaged.insert(corresp[nidx]);
         }
         auto nodes = topo->connectivity(dim0, dim1);
         // set damage field to MAX_DAM for nodes related to marked nodes
         // and count edges
         auto d_val = d->x()->mutable_array();
         p=0;
         for (auto &di : d_val)
         {
            if (node_damaged.find(p) != node_damaged.end())
            {
               di = MAX_DAM;
            }
            ++p;
         }
         // sum damage in owners
         d->x()->scatter_rev(std::plus<scalar>());
         // an owner may have damage above MAX_DAM if edges related to this owner are present in
         // more than one process: reset do MAX_DAM
         for (auto &di : d_val)
            if (di > MAX_DAM) di = MAX_DAM;
         // propagate to ghost: a ghost may still be null if no edges are marked in current process but in other process ghost or
         // owner counterparts are not null, thus this ghost has to be updated
         d->x()->scatter_fwd();

         // set edge counter
         int32_t nloc_edges = topo->index_map(dim1)->size_local();
         for (int32_t l = 0; l < nlg; ++l)
         {
            auto link = nodes->links(l);
            nd_val[l] = 0.;
            for (auto ll : link)
               // edge owned ? yes, count contribution
               if (ll < nloc_edges) ++nd_val[l];
         }
         nd_valv.scatter_rev(std::plus<scalar>());
         for (int32_t l = 0; l < nloc; ++l) num_edges_per_nodes[l] = 1./nd_val[l];

         // smooth
         std::int32_t id;
         for (short iter_smooth = 0, max_smooth = 8 * (MAX_REFINE + 1); iter_smooth < max_smooth; ++iter_smooth)
         {
            // first loop to enlarge
            for (int32_t l = 0; l < nlg; ++l)
            {
               auto &di = d_val[corresp[l]];
               if (di < 0.01)
               {
                  scalar dam = 0.;
                  auto link = nodes->links(l);
                  for (auto ll : link)
                  {
                     // only count contribution from owned edges to avoid duplication
                     if (ll < nloc_edges)
                     {
                        auto nn = edges->links(ll);
                        id = nn[0];
                        if (id == l) id = nn[1];
                        dam += d_val[corresp[id]];
                     }
                  }
                  nd_val[l] = dam;
               }
               else
                  nd_val[l] = 0.;
            }

            // sum all contribution on owner
            nd_valv.scatter_rev(std::plus<scalar>());

            // compute final value on owner into d_val (for owner dofs/mesh index are the same)
            for (int32_t l = 0; l < nloc; ++l) d_val[l] = std::max(nd_val[l] * num_edges_per_nodes[l], d_val[l]);

            // impose owner on ghost
            d->x()->scatter_fwd();

            // loop again to smooth and enlarge
            for (int32_t l = 0; l < nlg; ++l)
            {
               auto &di = d_val[corresp[l]];
               scalar dam = 0.;
               auto link = nodes->links(l);
               for (auto ll : link)
               {
                  // only count contribution from owned edges to avoid duplication
                  if (ll < nloc_edges)
                  {
                     auto nn = edges->links(ll);
                     id = nn[0];
                     if (id == l) id = nn[1];
                     dam += d_val[corresp[id]];
                  }
               }
               nd_val[l] = dam;
            }
            // sum, compute, impose
            nd_valv.scatter_rev(std::plus<scalar>());
            for (int32_t l = 0; l < nloc; ++l) d_val[l] = std::max(nd_val[l] * num_edges_per_nodes[l], d_val[l]);
            d->x()->scatter_fwd();
         }
      }
      measure[6]+=MPI_Wtime()-start; // 4.2 Define damage
      //d->x()->set(0.);
      // output damage to file 
      start = MPI_Wtime();
#ifdef USE_XDMF_FOR_OUTPUT
#ifdef USE_SQUARE
      std::string filename = "run/d_square_" + std::to_string(nb_proc);
      {
         auto v = d->x()->array();
         size_t i= 0;
         for (auto &di : v)
         cout<<i++<<" "<<di<<endl;
      }
#else
      std::string filename = "run/d_" + std::to_string(nb_proc);
#endif
      filename += ".xdmf";
      dolfinx::io::XDMFFile::Encoding  encoding = dolfinx::io::XDMFFile::Encoding::HDF5;
      //dolfinx::io::XDMFFile::Encoding  encoding = dolfinx::io::XDMFFile::Encoding::ASCII;
      {
         dolfinx::io::XDMFFile ofile(MPI_COMM_WORLD, filename, "w",encoding);
         ofile.write_mesh(*pmesh);
         ofile.write_function(*d, 0.);
         ofile.close();
      }
#endif

      measure[12]+=MPI_Wtime()-start; // 8 Outputs
      start = MPI_Wtime();


      //== create Young modulus field ================================
      auto E = std::make_shared<fem::Function<scalar>>(ES);
      auto part = std::make_shared<fem::Function<scalar>>(ES);

      //== specific field to output partition ========================
      part->name="partition";
      part->x()->set(mpi_rank);
#ifdef USE_XDMF_FOR_OUTPUT
      start = MPI_Wtime();
      // output partition
      filename = "run/partition_" + std::to_string(nb_proc) + ".xdmf";
      {
         dolfinx::io::XDMFFile ofile(MPI_COMM_WORLD, filename, "w");
         ofile.write_mesh(*pmesh);
         ofile.write_function(*part, 0.);
      }
      measure[12]+=MPI_Wtime()-start; // 8 Outputs
#endif

      //== fill Young modulus field ==================================
      start = MPI_Wtime();
      E->name="E";
      auto val_E=E->x()->mutable_array();
      // no orientation taking into account for now
      // all domains (i.e. unique gem) local to this proc are given by inspecting  gmsh physical tags
      // For each tag set an arbitrary Young modulus based on physical id
      // For now we consider 200 Young values from 5.e6 to 1.e8 with a semi-random variation
      std::vector<scalar> E_range(200);
      srand(6575); // fix seed so that E_range is the same across process and for all runs
      for (auto &v : E_range)
      {
         const scalar a=(1.e8-5.e6)/199.;
         v=a*(rand()%200)+5.e6;
#ifdef USE_ECST
         v = 1.e6;
#endif
      }
      auto local_tag = cell_tag.values();
      std::transform(local_tag.begin(), local_tag.end(), val_E.begin(),
                     [&E_range](int32_t phys) -> scalar { return E_range[phys % 200]; });
      measure[4]+=MPI_Wtime()-start; // 4.1 Material constant
#ifdef USE_XDMF_FOR_OUTPUT
      start = MPI_Wtime();
      // output E to file
      filename = "run/E_" + std::to_string(nb_proc) + ".xdmf";
      {
         dolfinx::io::XDMFFile ofile(MPI_COMM_WORLD, filename, "w");
         ofile.write_mesh(*pmesh);
         ofile.write_function(*E, 0.);
      }
      measure[12]+=MPI_Wtime()-start; // 8 Outputs
#endif

      //== volumetric load ===========================================
      start = MPI_Wtime();
      auto f = std::make_shared<fem::Function<scalar>>(V);
      f->name="f";
#ifdef USE_VOLUME
      f->interpolate([](auto x) -> std::pair<std::vector<scalar>, std::vector<std::size_t>> {
         size_t nb_col = x.extent(0);
         size_t nb_row = x.extent(1);
         size_t nb_col_v = 2;
         std::vector<scalar> vdata(nb_col_v * nb_row);
         namespace stdex = MDSPAN_IMPL_STANDARD_NAMESPACE::MDSPAN_IMPL_PROPOSED_NAMESPACE;
         MDSPAN_IMPL_STANDARD_NAMESPACE::mdspan<
             scalar, MDSPAN_IMPL_STANDARD_NAMESPACE::extents<std::size_t, 2, MDSPAN_IMPL_STANDARD_NAMESPACE::dynamic_extent>>
             v(vdata.data(), nb_col_v, nb_row);
         for (std::size_t p = 0; p < nb_row; ++p)
         {
            double xf = 100000.;
            //xf *= 1. / (0.0001 + x(0, p));
            double r=x(0, p)-0.5;
            xf *= -r * r * r;
            // double xf = x(0, p);
            double y = x(1, p) - 0.5;
            v(0, p) = (1600. * y * y - 500.) * xf;
            v(1, p) = 0.;
         }
         return {vdata, {nb_col_v, nb_row}};
      });
#else
      f->x()->set(0.);
#endif
      //== surface load marking ======================================
      std::map<dolfinx::fem::IntegralType, std::vector<std::pair<std::int32_t, std::span<const std::int32_t>>>> subdomains
#ifdef USE_SURF
      ;
      const std::vector<std::int32_t> upper_lower = dolfinx::mesh::locate_entities_boundary(*pmesh, dim1, [&eps,&eps1](auto x) {
         std::vector<int8_t> marker(x.extent(1), false);
         for (std::size_t p = 0; p < x.extent(1); ++p)
         {
            double y = x(1, p);
            if (y < eps || y>eps1) marker[p] = true;
         }
         return marker;
      });
      subdomains.insert(
          std::make_pair(dolfinx::fem::IntegralType::exterior_facet,
                         std::vector<std::pair<std::int32_t, std::span<const std::int32_t>>>(
                             1, std::make_pair(0, std::span<const std::int32_t>(upper_lower.data(), upper_lower.size())))));
#else
      ={};
#endif

          measure[8] += MPI_Wtime() - start;  // 5.2 Neuman setting

      //== create function to hold solution ==========================
      auto u = std::make_shared<dolfinx::fem::Function<scalar>>(V);
      string n_disp = "disp_FeniCSx_";
      n_disp += SEXT;
      u->name = n_disp;
      // innit to zero
      u->x()->set(0.);

      //== imposed boundary condition : Dirichlet ====================
      start = MPI_Wtime();
      // left edge fixed to null value in all direction (x,y and z)
      auto u0 = std::make_shared<dolfinx::fem::Function<scalar>>(V);
      u0->x()->set(0);
      // Find nodes with bc applied (order=1 => only nodes count => in // safe, no need to communicate !!ERROR!! By rereading doc it is
      // not so clear. If locate_entities_boundary use face/edge to select nodes its wrong)
      const std::vector<std::int32_t> bc_nodes = dolfinx::mesh::locate_entities_boundary(*pmesh, dim0, [&eps](auto x) {
         std::vector<int8_t> marker(x.extent(1), false);
         for (std::size_t p = 0; p < x.extent(1); ++p)
         {
            double x0 = x(0, p);
            if (std::abs(x0) < eps) marker[p] = true;
         }
         return marker;
      });
      // Find constrained dofs
      const std::vector<std::int32_t> bdofs =
          dolfinx::fem::locate_dofs_topological(*V->mesh()->topology_mutable(), *V->dofmap(), dim0, bc_nodes);
      // Fixed nodes bc
      auto bcl = std::make_shared<dolfinx::fem::DirichletBC<scalar>>(u0, bdofs);
      // right edge fixed to x=-+0.01 value 
#ifdef USE_TRAC
      scalar imp[3] = {0.01, 0.};
#else
      scalar imp[3] = {-0.01, 0.};
#endif
      const dolfinx::fem::Constant<scalar> imp_cst(std::span<const scalar>(imp, 2));
      auto impdisp =
          std::make_shared<const dolfinx::fem::Constant<scalar>>(dolfinx::fem::Constant<scalar>(std::span<const scalar>(imp, 2)));
      // Find nodes with bc applied (order=1 => only nodes count => in // safe, no need to communicate)
      const std::vector<std::int32_t> imp_nodes = dolfinx::mesh::locate_entities_boundary(*pmesh, dim0, [&eps](auto x) {
         std::vector<int8_t> marker(x.extent(1), false);
         for (std::size_t p = 0; p < x.extent(1); ++p)
         {
            double x0 = x(0, p);
            if (std::abs(x0-1.) < eps) marker[p] = true;
         }
         return marker;
      });
      // Find constrained dofs
      const std::vector<std::int32_t> idofs =
          dolfinx::fem::locate_dofs_topological(*V->mesh()->topology_mutable(), *V->dofmap(), dim0, imp_nodes);
      // impose nodes bc
      auto bcr = std::make_shared<dolfinx::fem::DirichletBC<scalar>>(impdisp, idofs, V);
      // init u with bc via bcr (bcl gives only zero so not used here)
      // With this setting non-linear loop start at least with non-zero imposed displacement value set
      // bcr->set(u->x()->mutable_array());
      // Does not have impact on performance thus it is not used
      measure[7]+=MPI_Wtime()-start; // 5.2 Dirichlet setting
      start = MPI_Wtime();

      //== Coefficient field associated to forms =====================
      // UFC E (Young modulus) is defined in  python form
      const std::map<std::string, std::shared_ptr<const dolfinx::fem::Function<scalar>>> coefficients = {
          {"E", E}, {"d", d}, {"f", f}, {"u", u}};

      //== Form to compute residual ==================================
      // UFC FORM_F comes from c  form file generated from python with ffcx
      auto F_form = std::make_shared<dolfinx::fem::Form<scalar, double>>(
          dolfinx::fem::create_form<scalar>(*FORM_F, {V}, coefficients, constants, subdomains));

      //== Form to compute residual ==================================
      // UFC FORM_J comes from c form file generated from python with ffcx
      auto J_form = std::make_shared<dolfinx::fem::Form<scalar, double>>(
          dolfinx::fem::create_form<scalar>(*FORM_J, {V,V}, coefficients, constants, {}));

      //== non Linear matrix related to Jacobian form ================
      auto A = std::make_shared<dolfinx::la::petsc::Matrix>(dolfinx::fem::petsc::create_matrix(*J_form), false);

      //== residual vector ===========================================
      dolfinx::la::Vector<scalar> b(V->dofmap()->index_map,
                                    V->dofmap()->index_map_bs());
      b.set(0);

      //== Warp vector in petsc vector ===============================
      dolfinx::la::petsc::Vector b_petsc(dolfinx::la::petsc::create_vector_wrap(b), false);
      dolfinx::la::petsc::Vector u_petsc(dolfinx::la::petsc::create_vector_wrap(*(u->x())), false);

      //== Warp u in a span ==========================================
      auto u_array=u->x()->mutable_array();
      measure[9]+=MPI_Wtime()-start; // 7.1 Nonlinear form creation
      start = MPI_Wtime();

      //== non linear solver creation ================================
      dolfinx::nls::petsc::NewtonSolver newton_solver(MPI_COMM_WORLD);

      //== non linear solver setting =================================
      // max iter
      newton_solver.max_it=10;
      newton_solver.error_on_nonconvergence=true;
      // tolerance
      newton_solver.rtol = 1.e-7;
      newton_solver.atol = 5.e-8;
      newton_solver.report=true;

      // setting ksp tolerance
      KSP ksp = newton_solver.get_krylov_solver().ksp();
      KSPSetTolerances(ksp, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT, 2000);

#ifdef USE_BOOMERAMG
      // use cg and hypre boomeramg preconditioner to compare to mfem (can't be changed by user) 
      PetscOptionsSetValue(NULL, "-nls_solve_ksp_type", "cg");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_type", "hypre");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_type", "boomeramg");
      //PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_print_debug", "1");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_print_statistics", "1");
#ifdef USE_BOOMERAMG_FINE_TUNNING
#if (USE_BOOMERAMG_FINE_TUNNING == 0)
      ////// OLD ATTEMPTED //////////////////////////////
      ////// Try to mimic SetElasticityOptions but some tests with MFEM show that it is
      ////// not the best.
      ////// This is an attempt has it is unstable with many proc and non null MAXREFINE 
      ///////////////////////////////////////////////////
      // fine tuning BoomerAMG hypre precond (can't be changed by user)
      // comes in particular from MFEM SetElasticityOptions
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_coarsen_type", "HMIS");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_nodal_coarsen", "4");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_nodal_coarsen_diag", "1");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_relax_type_up", "l1scaled-SOR/Jacobi");   // "8"
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_relax_type_down", "l1scaled-SOR/Jacobi"); // "8"
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_relax_type_coarse", "l1scaled-SOR/Jacobi"); // "8"
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_vec_interp_variant", "2");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_vec_interp_qmax", "4");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_vec_interp_smooth", "1");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_interp_refine", "1");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_interp_type", "ext+i");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_strong_threshold", "0.5");
      // reset all from setting above or user option
      KSPSetFromOptions(ksp);
      // Near null space for boomeramg
      // Comes from MFEM again (RecomputeRBMs)
      // Here in 2D only one vector is created. It corresponds to rigid body rotation (translation in x and y are not set)
      auto rbm = std::make_shared<fem::Function<scalar>>(V);
      rbm->interpolate([](auto x) -> std::pair<std::vector<scalar>, std::vector<std::size_t>> {
         size_t nb_col = x.extent(0);
         size_t nb_row = x.extent(1);
         size_t nb_col_v = 2;
         std::vector<scalar> vdata(nb_col_v * nb_row);
         namespace stdex = MDSPAN_IMPL_STANDARD_NAMESPACE::MDSPAN_IMPL_PROPOSED_NAMESPACE;
         MDSPAN_IMPL_STANDARD_NAMESPACE::mdspan<
             scalar, MDSPAN_IMPL_STANDARD_NAMESPACE::extents<std::size_t, 2, MDSPAN_IMPL_STANDARD_NAMESPACE::dynamic_extent>>
             v(vdata.data(), nb_col_v, nb_row);
         double c = 1. / sqrt(2.);  // scaling factor
         for (std::size_t p = 0; p < nb_row; ++p)
         {
            double xv = x(0, p) - 0.5;
            double yv = x(1, p) - 0.5;
            v(0, p) = (xv + yv) * c - xv;
            v(1, p) = (yv - xv) * c - yv;
         }
         return {vdata, {nb_col_v, nb_row}};
      });
#if 0
      {
         rbm->name = "rbm";
         std::string filename = "run/rbm_" + std::to_string(nb_proc);
         filename += ".xdmf";
         dolfinx::io::XDMFFile::Encoding encoding = dolfinx::io::XDMFFile::Encoding::HDF5;
         dolfinx::io::XDMFFile ofile(MPI_COMM_WORLD, filename, "w", encoding);
         ofile.write_mesh(*pmesh);
         ofile.write_function(*rbm, 0.);
         ofile.close();
      }
#endif
      dolfinx::la::petsc::Vector rbm_petsc(dolfinx::la::petsc::create_vector_wrap(*(rbm->x())), false);
      const Vec basis = rbm_petsc.vec();
      // create matNullSpace object associated to basis
      MatNullSpace matNullSpace;
      MatNullSpaceCreate(MPI_COMM_WORLD, PETSC_FALSE, 1, &basis, &matNullSpace);
      // MatNullSpaceCreate(MPI_COMM_WORLD, PETSC_TRUE, 1, &basis, &matNullSpace);

      // add matNullSpace object to Mat A matrix so that during hypre boomeramg pc creation nullspace basis is
      // taken into account
      // Attached to A and as NewtonSolver::setP is not use, A is also used for PC. Thus  matNullSpace should be discover during
      // PCSetup
      Mat Amat = A->mat();
      MatSetNearNullSpace(Amat, matNullSpace);
      MatNullSpaceDestroy(&matNullSpace);
#else
      ////// FOLLOW MANUAL MFEM TUNNING /////////////////
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_coarsen_type", "HMIS");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_relax_type_up", "l1scaled-SOR/Jacobi");   // "8"
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_relax_type_down", "l1scaled-SOR/Jacobi"); // "8"
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_relax_type_coarse", "Gaussian-elimination"); // "9"
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_interp_type", "ext+i");
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_numfunctions", "2"); // Number of functions set to 2 as in 2D
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_agg_nl", "0"); // No aggressive coarsening
      PetscOptionsSetValue(NULL, "-nls_solve_pc_hypre_boomeramg_strong_threshold", "0.25");
#endif
#endif
#endif
      // reset all from user option with above setting masking user's eventual choice
      KSPSetFromOptions(ksp);

      measure[10]+=MPI_Wtime()-start; // 7.2 Solver creation
      // Function to compute residual vector
      newton_solver.setF(
          [&F_form, &J_form, &b, &bcl, &bcr, &u_array,&measure](const Vec x_petsc, Vec b_petsc) {
             double s=MPI_Wtime();
             // Directly use b captured instead of petsc b_petsc argument to use dolfinx interface
             // reset to null
             b.set(0.0);

             // assemble residual form
             dolfinx::fem::assemble_vector<scalar>(b.mutable_array(), *F_form);
             // compute Dirichlet contribution to the residual vector via J form (-k_if.u_i)
             dolfinx::fem::apply_lifting<scalar, double>(b.mutable_array(), {J_form}, {{bcl,bcr}}, {u_array}, -1.);

             // simpler to call b.scatter_rev(std::plus<scalar>()); but do not reset Vec petsc norm !
             VecGhostUpdateBegin(b_petsc, ADD_VALUES, SCATTER_REVERSE);
             VecGhostUpdateEnd(b_petsc, ADD_VALUES, SCATTER_REVERSE);

             // assign fixed values variation to rhs b
             // As Dirichlet dof in matrix J is an identity matrix rhs value becomes solution values and thus corresponds to
             // imposed fixed values variation
             dolfinx::fem::set_bc<scalar, double>(b.mutable_array(), {bcl, bcr},u_array,-1.);

             // force usage of b_petsc so that norm get update
             // add 0 to arbitrary first component
             //scalar* array = nullptr;
             //VecGetArray(b_petsc, &array);
             //VecRestoreArray(b_petsc, &array);
             measure[13]+=MPI_Wtime()-s; // 6.3 Create and assemble elementary vector
          },
          b_petsc.vec());
      // Function to compute Jacobian matrix
      newton_solver.setJ(
          [&J_form, &bcl, &bcr, &V, &measure](const Vec x_unused, Mat A_) {
             double s=MPI_Wtime();
             MatZeroEntries(A_);
             // all terms but with null row for Dirichlet dofs
             dolfinx::fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A_, ADD_VALUES), *J_form, {bcl, bcr});
             MatAssemblyBegin(A_, MAT_FLUSH_ASSEMBLY);
             MatAssemblyEnd(A_, MAT_FLUSH_ASSEMBLY);

             // set diagonal terms to 1 for Dirichlet dofs
             dolfinx::fem::set_diagonal<scalar>(la::petsc::Matrix::set_fn(A_, INSERT_VALUES), *V, {bcl, bcr}, 1.);
             MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
             MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
             measure[14]+=MPI_Wtime()-s; // 6.3 Create and assemble elementary matrix
          },
          A->mat());
      // Function to update in parrallel u field
      newton_solver.set_form([](Vec u_petsc) {
         VecGhostUpdateBegin(u_petsc, INSERT_VALUES, SCATTER_FORWARD);
         VecGhostUpdateEnd(u_petsc, INSERT_VALUES, SCATTER_FORWARD);
      });
      // Function to check convergence (verbose version of the default function)
      newton_solver.set_convergence_check(
            [&do_print](const nls::petsc::NewtonSolver& solver,
               const Vec r)->std::pair<double, bool>
            {
            scalar residual = 0.0;
            VecNorm(r, NORM_2, &residual);

            // Relative residual
            const double relative_residual = residual / solver.residual0();

            // Output iteration number and residual
            if (do_print)
               cout << scientific << "Newton iteration " << solver.iteration()
               << ": r (abs) = " << residual << " (tol = " << solver.atol
               << ") r (rel) = " << relative_residual << "(tol = " << solver.rtol
               << ")"<<endl;

            // Return true if convergence criterion is met
            if (relative_residual < solver.rtol or residual < solver.atol)
               return {residual, true};
            else
               return {residual, false};
            });
      //== non linear resolution =====================================
      start = MPI_Wtime();
      std::pair<int, bool> r = newton_solver.solve(u_petsc.vec());
      measure[11]+=MPI_Wtime()-start; // 7.3 NonLinear resolution
      KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
      if ( r.second )
      {
         if(do_print)
            cout<<"NL converge in "<<r.first<<" iteration(s)"<<endl;
      }
      else
      {
         if(do_print)
            cout<<"NL did not converge after "<<r.first<<" iteration(s)"<<endl;
         MPI_Abort(MPI_COMM_WORLD,-1234);
      }

      start = MPI_Wtime();
      //== strain field ==============================================
      constexpr auto family = basix::element::family::P;
#ifdef FORCE_NEW
      auto cell_types = pmesh->topology()->cell_types();
      auto cell_type = mesh::cell_type_to_basix_type(cell_types[0]);
#else
      auto cell_type = mesh::cell_type_to_basix_type(pmesh->topology()->cell_type());
#endif
      constexpr int k = 0;
      constexpr bool discontinuous = true;
      basix::FiniteElement S_element = basix::create_element<scalar_dolf>(family, cell_type, k, basix::element::lagrange_variant::unset,
                                                                basix::element::dpc_variant::unset, discontinuous);
      auto S = std::make_shared<fem::FunctionSpace<scalar_dolf>>(
          fem::create_functionspace(pmesh, S_element, std::vector<std::size_t>{3}));
      auto strain_expression = fem::create_expression<scalar, scalar_dolf>(*EXPR0_F, {{"u", u}}, {});

      auto strain = std::make_shared<fem::Function<scalar>>(S);
      string n_stra = "strain_FeniCSx_";
      n_stra += SEXT;
      strain->name = n_stra;
      strain->interpolate(strain_expression,*pmesh);
      std::string filename_strain;

      //== stress field ==============================================
      auto stress_expression = fem::create_expression<scalar, scalar_dolf>(*EXPR1_F, coefficients, constants);

      auto stress = std::make_shared<fem::Function<scalar>>(S);
      string n_stre = "stress_FeniCSx_";
      n_stre += SEXT;
      stress->name = n_stre;
      stress->interpolate(stress_expression,*pmesh);
      std::string filename_stress;
      measure[15]+=MPI_Wtime()-start; // 8.1 strain/stress computation
      start = MPI_Wtime();


      //== output solution to file ===================================

#ifdef USE_XDMF_FOR_OUTPUT
#ifdef USE_SQUARE
      filename = "run/square_damage_" + std::to_string(nb_proc) + ".xdmf";
      filename_strain = "run/square_damage_strain_" + std::to_string(nb_proc) + ".xdmf";
      filename_stress = "run/square_damage_stress_" + std::to_string(nb_proc) + ".xdmf";
#else
      string ext = SEXT;
      ext += "_";
#ifdef USE_TRAC
#ifdef USE_VOLUME
      filename = "run/asym_elasto_damage_trac_vol_" + ext + std::to_string(nb_proc) + ".xdmf";
      filename_strain = "run/asym_elasto_damage_strain_trac_vol_" + ext + std::to_string(nb_proc) + ".xdmf";
      filename_stress = "run/asym_elasto_damage_stress_trac_vol_" + ext + std::to_string(nb_proc) + ".xdmf";
#else
      filename = "run/asym_elasto_damage_trac_" + ext + std::to_string(nb_proc) + ".xdmf";
      filename_strain = "run/asym_elasto_damage_strain_trac_" + ext + std::to_string(nb_proc) + ".xdmf";
      filename_stress = "run/asym_elasto_damage_stress_trac_" + ext + std::to_string(nb_proc) + ".xdmf";
#endif
#else
      filename = "run/asym_elasto_damage_comp_" + ext + std::to_string(nb_proc) + ".xdmf";
      filename_strain = "run/asym_elasto_damage_strain_comp_" + ext + std::to_string(nb_proc) + ".xdmf";
      filename_stress = "run/asym_elasto_damage_stress_comp_" + ext + std::to_string(nb_proc) + ".xdmf";
#endif
#endif
      dolfinx::io::XDMFFile ofile(MPI_COMM_WORLD, filename, "w");
      ofile.write_mesh(*pmesh);
      ofile.write_function(*u, 0.0);

      dolfinx::io::XDMFFile ofile_strain(MPI_COMM_WORLD, filename_strain, "w");
      ofile_strain.write_mesh(*pmesh);
      ofile_strain.write_function(*strain, 0.0);

      dolfinx::io::XDMFFile ofile_stress(MPI_COMM_WORLD, filename_stress, "w");
      ofile_stress.write_mesh(*pmesh);
      ofile_stress.write_function(*stress, 0.0);
#elif defined(USE_VTK_FOR_OUTPUT)
      dolfinx::io::VTKFile ofile(MPI_COMM_WORLD, filename, "w");
      ofile.write(*pmesh);
      // to be finished
      std::vector<std::reference_wrapper<const fem::Function<scalar, scalar_dolf>>> func={reference_wrapper(*u)};
      ofile.write(func, 0.0);
      ofile.close();
#elif defined(USE_ADIOS_FOR_OUTPUT)
#ifdef USE_SQUARE
      std::string filename = "run/square_damage_";
#else
      string ext = SEXT;
      ext += "_";
#ifdef USE_TRAC
#ifdef USE_VOLUME
      std::string filename = "run/asym_elasto_damage_trac_vol_" + ext;
#else
      std::string filename = "run/asym_elasto_damage_trac_" + ext;
#endif
#else
      std::string filename = "run/asym_elasto_damage_comp_" + ext;
#endif
#endif
      std::string f1 = filename + "vect_" + std::to_string(nb_proc) + ".bp";
      {
         dolfinx::io::adios2_writer::U<scalar_dolf> fields = {u, f};
         dolfinx::io::VTXWriter<scalar_dolf> vtx(MPI_COMM_WORLD, f1, fields, "bp4");
         vtx.write(0);
      }
      std::string f2 = filename + "DGscal_" + std::to_string(nb_proc) + ".bp";
      {
         dolfinx::io::adios2_writer::U<scalar_dolf> fields = {E, part};
         dolfinx::io::VTXWriter<scalar_dolf> vtx(MPI_COMM_WORLD, f2, fields, "bp4");
         vtx.write(0);
      }
      std::string f3 = filename + "scal_" + std::to_string(nb_proc) + ".bp";
      {
         dolfinx::io::adios2_writer::U<scalar_dolf> fields = {d};
         dolfinx::io::VTXWriter<scalar_dolf> vtx(MPI_COMM_WORLD, f3, fields, "bp4");
         vtx.write(0);
      }
#if 1
      std::string f4 = filename + "tens_" + std::to_string(nb_proc) + ".bp";
      {
         dolfinx::io::adios2_writer::U<scalar_dolf> fields = {strain, stress};
         dolfinx::io::VTXWriter<scalar_dolf> vtx(MPI_COMM_WORLD, f4, fields, "bp4");
         vtx.write(0);
      }
#endif
#endif
      measure[12]+=MPI_Wtime()-start; // 8 Outputs
      measure[0] += MPI_Wtime() - start_all;  // All


#ifdef IN_COMP
      //== field comparison with mfem ================================
      { //change scope to isolate this part
         assert(nb_proc == 1);
#ifdef USE_VOLUME
         string in_file_name = "data/mfem_disp_"+std::to_string(MAX_REFINE)+"_vol";
         std::string f1 = filename + "diff_" + std::to_string(MAX_REFINE) + "_vol.bp";
         std::string f2 = filename + "energ_" + std::to_string(MAX_REFINE) + "_vol.bp";
#else
         string in_file_name = "data/mfem_disp_"+std::to_string(MAX_REFINE);
         std::string f1 = filename + "diff_" + std::to_string(MAX_REFINE) + ".bp";
         std::string f2 = filename + "energ_" + std::to_string(MAX_REFINE) + ".bp";
#endif
         ifstream infile(in_file_name,std::ifstream::binary);
         double cx, cy, dx, dy;
         size_t siz=sizeof(double);
         // read raw data
         container_t disp_mfem;
         infile.read(reinterpret_cast<char *>(&cx), siz);
         while (infile.good())
         {
            infile.read(reinterpret_cast<char *>(&cy), siz);
            infile.read(reinterpret_cast<char *>(&dx), siz);
            infile.read(reinterpret_cast<char *>(&dy), siz);
            disp_mfem.push_back({cx, cy, dx, dy});
            infile.read(reinterpret_cast<char *>(&cx), siz);
         }
         //sort by cx,cy
         double eps=EPS;
         short ind;
         auto lti = [&eps, &ind](const quarted_t &left, const quarted_t &right) -> bool {
            double diff = left[ind] - right[ind];
            //cout<<" "<<ind<<" "<<left[ind]<<" "<<right[ind]<<" "<<diff;
            if (!double_eq(left[ind], 0.) && !double_eq(right[ind], 0.)) diff /= fabs(left[ind]);
            //cout<<" "<<diff<<" "<<(diff < -eps)<<endl;
            if (diff < -eps) return true;
            return false;
         };
         auto lt=[&ind,&lti](const quarted_t &left, const quarted_t &right) -> bool {
            ind=0;
            if (lti(left,right)) return true;
            if (lti(right,left)) return false;
            ind=1;
            return lti(left,right);
         };
         std::sort(disp_mfem.begin(), disp_mfem.end(), lt);
         //for (auto &v : disp_mfem) cout << "mfem " << v[0] << " " << v[1] << endl;

         auto tdc = V->tabulate_dof_coordinates(false);
         size_t stdc=tdc.size();
         eps=1.e-5;
         double errx=0.;
         double erry=0.;
         for (size_t p = 0, q = 0; p < stdc; p += 3)
         {
            cx = tdc[p];
            cy = tdc[p + 1];
            quarted_t cur = {cx, cy, 0., 0.};
            ind = 0;
            auto itbx = lower_bound(disp_mfem.begin(), disp_mfem.end(), cur, lti);
            bool notfound=true;
            dx=dy=0.;
            while (itbx != disp_mfem.end())
            {
               auto nd = *itbx;
               //cout << " " << nd[0] << " " << nd[1] << " " << double_eq(cx, nd[0], eps) << " " << double_eq(cy, nd[1], eps) << endl;
               if (double_eq(cx, nd[0], eps))
               {
                  if (double_eq(cy, nd[1], eps))
                  {
                     dx = nd[2];
                     dy = nd[3];
                     notfound = false;
                     break;
                  }
               }
               else
                  break;
               ++itbx;
            }
            if (notfound)
            {
               cout << "dof " << p / 3;
               cout << " " << cx << " " << cy;
               cout << " notfound"<<endl;
               throw -123;
            }
            u_array[q] -= dx;
            u_array[q + 1] -= dy;
            // compute norme per component
            errx+=u_array[q]*u_array[q];
            erry+=u_array[q+1]*u_array[q+1];
            q += 2;
         }
         double L2x=sqrt(errx);
         double L2y=sqrt(erry);
         cout<<"Error L2 x:"<<L2x<<endl;
         cout<<"Error L2 y:"<<L2y<<endl;

         auto energy_expression = fem::create_expression<scalar, scalar_dolf>(*EXPR2_F, coefficients, constants);
         auto energy_error = std::make_shared<fem::Function<scalar>>(ES);
         string n_ener = "Relative energy error strain_FeniCSx_";
         n_ener += SEXT;
         energy_error->name = n_ener;
         energy_error->interpolate(energy_expression, *pmesh);
         auto energy_error_array = energy_error->x()->mutable_array();
         double energy_error_sum = std::accumulate(energy_error_array.begin(), energy_error_array.end(), 0.);
         cout<<"Error in energy:"<<energy_error_sum<<endl;
         energy_error_sum = 100. / energy_error_sum;
         std::transform(energy_error_array.begin(),energy_error_array.end(),energy_error_array.begin(),[&energy_error_sum](const double &v){return v*energy_error_sum;});
         L2x = 100. / L2x;
         L2y = 100. / L2y;
         int q=0;
         std::transform(u_array.begin(), u_array.end(), u_array.begin(), [&L2x, &L2y, &q](const double &v) {
            if ((q++) % 2) return v * L2y;
            return v * L2x;
         });
#ifdef USE_ADIOS_FOR_OUTPUT
         {
            {
               u->name = "FEniCSx-MFEM";
               dolfinx::io::adios2_writer::U<scalar_dolf> fields = {u};
               dolfinx::io::VTXWriter<scalar_dolf> vtx(MPI_COMM_WORLD, f1, fields, "bp4");
               vtx.write(0);
            }
            {
               dolfinx::io::adios2_writer::U<scalar_dolf> fields = {energy_error};
               dolfinx::io::VTXWriter<scalar_dolf> vtx(MPI_COMM_WORLD, f2, fields, "bp4");
               vtx.write(0);
            }
         }
#endif
      }
#endif
      start_all = MPI_Wtime();
   }

   measure[0] += MPI_Wtime() - start_all;  // All
#ifdef PROF_KERNEL
   measure[16] = gmeasure1;  // kernel vect
   measure[17] = gmeasure2;  // kernel matrix
#endif
   double pow_measure[NB_MEASURE];
   double ecart_type_measure[NB_MEASURE];
   double max_measure[NB_MEASURE];
   double min_measure[NB_MEASURE];
   double avg_measure[NB_MEASURE];
   std::transform(measure,measure+NB_MEASURE,pow_measure,[](const double &v){return v*v;});
   MPI_Reduce(measure, max_measure, NB_MEASURE, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   MPI_Reduce(measure, min_measure, NB_MEASURE, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
   MPI_Reduce(measure, avg_measure, NB_MEASURE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(pow_measure, ecart_type_measure, NB_MEASURE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if (do_print)
   {
      cout << setfill('=') << setw(137) << " " << endl;
      cout << setfill(' ') << "| " << setw(12) << "val min"
           << " | " << setw(12) << "val max"
           << " | " << setw(12) << "std dev"
           << " | " << setw(12) << "% Coef var"
           << " | " << setw(12) << "val avg"
           << " | " << setw(12) << "% total avg"
           << " | " << setw(44) << " |" << endl;
      SL(0, "All")
      SL(1, "1. Initialize ")
      SL(2, "2.1 Read the mesh")
      SL(3, "2.2 Refine the mesh")
      SL(5, "3.1  Define space")
      SL(6, "3.2  Define damage")
      SL(4, "4.1  Material constant")
      SL(7, "5.1 Dirichlet setting")
      SL(8, "5.2 Neuman setting");
      SL(13, "6.3 Create and assemble elementary vector");
      SL(14, "6.4 Create and assemble elementary matrix");
#ifdef PROF_KERNEL
      SL(16, "6.3 Create Elementary vector");
      SL(17, "6.4 Create Elementary matrix");
#endif
      SL(9, "7.1 Nonlinear form creation");
      SL(10, "7.2 Solver creation");
      SL(11, "7.3 NonLinear resolution");
      SL(12, "8 Outputs");
      SL(15, "8.1 strain/stress computation");
      cout << setfill('=') << setw(137) << " " << endl;
   }

   //== clean/close/exit ======================================
   // close petsc environment
   // imply close distributed environment
   PetscFinalize();
   return 0;
}
