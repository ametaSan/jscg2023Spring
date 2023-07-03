////////////////////////////////////////////////////////////////////
//
// $Id: LoopSubL.hxx 2022/05/23 01:00:35 kanai Exp $
//
// Copyright (c) 2022 Takashi Kanai
// Released under the MIT license
//
////////////////////////////////////////////////////////////////////

#ifndef _LOOPSUB_HXX
#define _LOOPSUB_HXX 1

#include "envDep.h"
#include "mydef.h"

#include <cmath>
#include <vector>
using namespace std;

#include "myEigen.hxx"

#include "EdgeL.hxx"
#include "FaceL.hxx"
#include "MeshL.hxx"
#include "VertexL.hxx"
#include "VertexLCirculator.hxx"
#include "MeshUtiL.hxx"

// masks
#define LOOP_MASK_BETA3 0.1875
#define LOOP_MASK_12 0.5
#define LOOP_MASK_18 0.125
#define LOOP_MASK_38 0.375
#define LOOP_MASK_58 0.625
#define LOOP_MASK_68 0.75

// Loop beta for valence = 4 ... 10
static double loop_beta[7] = {
  0.121094,
  0.0840932,
  0.0625,
  0.0490249,
  0.0400678,
  0.033785,
  0.0291778
};


class LoopSub {
 public:
  LoopSub() : mesh_(NULL), submesh_(NULL){};
  LoopSub(MeshL& mesh) : submesh_(NULL) { setMesh(mesh); };
  LoopSub(MeshL& mesh, MeshL& submesh) {
    setMesh(mesh);
    setSubMesh(submesh);
  };
  ~LoopSub(){};

  void setMesh(MeshL& mesh) { mesh_ = &mesh; };
  void setSubMesh(MeshL& mesh) { submesh_ = &mesh; };
  MeshL& submesh() const { return *submesh_; };

  bool emptyMesh() const { return (mesh_ != NULL) ? false : true; };
  bool emptySubMesh() const { return (submesh_ != NULL) ? false : true; };

  void apply() {
    if (init() == false) return;
    setSplit();
    setStencil();
    submesh_->calcAllFaceNormals();
    std::cout << "loop subdiv.: done. v " << submesh_->vertices_size()
              << " f " << submesh_->faces_size() << std::endl;
    clear();
  };

  bool init() {
    if (emptyMesh()) return false;
    if (emptySubMesh()) return false;

    // check for rectangle faces
    for (auto fc : mesh_->faces()) {
      if (fc->size() != TRIANGLE) {
        std::cerr << "Error: A non-triangle face is included. " << std::endl;
        return false;
      }
    }

    // don't delete edges
    mesh_->createConnectivity(false);

    return true;
  };

  void clear() {
    even_.clear();
    odd_.clear();
  };

  // ここに分割処理のコードを追加してください．
  void setSplit() 
  {
    Eigen::Vector3d p = Eigen::Vector3d::Zero();

    // submesh の even vertex を生成 (mesh の頂点数)
    even_.resize(mesh_->vertices_size());
    std::cout << "vertices size: " << mesh_->vertices_size() << std::endl;
    for (int i = 0; i < mesh_->vertices_size(); ++i)
        even_[i] = submesh_->addVertex(p);

    // submesh の odd vertex を作成 (mesh のエッジの数)
    odd_.resize(mesh_->edges().size());
    std::cout << "edges size: " << mesh_->edges().size() << std::endl;
    for (int i = 0; i < mesh_->edges().size(); ++i)
        odd_[i] = submesh_->addVertex(p);

    for (auto fc : mesh_->faces()) 
    {
        // ここに 4-to-1 分割の新しい4つの面を submesh に作成するコードを追加
        HalfedgeL* he = fc->begin();
        VertexL* v1 = even_[he->vertex()->id()];
        VertexL* v2 = even_[he->next()->vertex()->id()];
        VertexL* v3 = even_[he->prev()->vertex()->id()];
        VertexL* e1 = odd_[he->edge()->id()];
        VertexL* e2 = odd_[he->next()->edge()->id()];
        VertexL* e3 = odd_[he->prev()->edge()->id()];

        FaceL* fc0 = submesh_->addFace();
        submesh_->addHalfedge(fc0, e1);
        submesh_->addHalfedge(fc0, e2);
        submesh_->addHalfedge(fc0, e3);

        FaceL* fc1 = submesh_->addFace();
        submesh_->addHalfedge(fc1, v1);
        submesh_->addHalfedge(fc1, e1);
        submesh_->addHalfedge(fc1, e3);

        FaceL* fc2 = submesh_->addFace();
        submesh_->addHalfedge(fc2, v2);
        submesh_->addHalfedge(fc2, e2);
        submesh_->addHalfedge(fc2, e1);

        FaceL* fc3 = submesh_->addFace();
        submesh_->addHalfedge(fc3, v3);
        submesh_->addHalfedge(fc3, e3);
        submesh_->addHalfedge(fc3, e2);
    }
  };


  // ここに位置計算処理のコードを追加してください．
  void setStencil() 
  {
    // even vertex
    for (auto vt : mesh_->vertices()) 
    {
      Eigen::Vector3d p;
      Eigen::Vector3d self_posi = vt->point();
      
      // 頂点に隣接する頂点の探索
      VertexLCirculator vc(vt);
      VertexL* vv = vc.beginVertexL();
      int ver_counter = 0;  // counting #vertices
      Eigen::Vector3d vertices_posi[20]; // 20 may be a good limit of #vertices
      do
      {
          vertices_posi[ver_counter] = vv->point(); //i don't like use referenced value
          ver_counter++;
          vv = vc.nextVertexL();
      } while ((vv != vc.firstVertexL()) && (vv != NULL));

      // loop subdivision
      if (ver_counter < 3)
      {
          std::cout << "Problems occur in vertices counting." << std::endl;
      }
      else if (ver_counter == 3)
      {
          Eigen::Vector3d sum_posi(0, 0, 0);
          for (int i = 0; i < ver_counter; i++)
          {
              sum_posi += vertices_posi[i];
          }
          p = self_posi * 0.4375f + sum_posi * 0.1875f;
      }
      else
      {
          Eigen::Vector3d sum_posi(0, 0, 0);
          for (int i = 0; i < ver_counter; i++)
          {
              sum_posi += vertices_posi[i];
          }
          float root_temp = 0.375f + (0.25f * cos(2.f * M_PI / ver_counter));
          float w_nextPosi = (1.f / ver_counter) * (0.625f - root_temp * root_temp);

          p = self_posi * (1.f - (ver_counter * w_nextPosi)) + sum_posi * w_nextPosi;
      }

      even_[vt->id()]->setPoint(p);
    }

    // odd vertex points
    for (auto ed : mesh_->edges()) 
    {
      Eigen::Vector3d p;
      // ここに odd vertex の頂点位置の計算コードを追加
      Eigen::Vector3d v1 = ed->lhe()->vertex()->point();
      Eigen::Vector3d v2 = ed->rhe()->vertex()->point();
      Eigen::Vector3d v3 = ed->lhe()->prev()->vertex()->point();
      Eigen::Vector3d v4 = ed->rhe()->prev()->vertex()->point();
      p = 0.375f * v1 + 0.375f * v2 + 0.125f * v3 + 0.125f * v4;

      odd_[ed->id()]->setPoint(p);
    }
  };

  double beta(int valence) {
    if (valence == 3) 
      return LOOP_MASK_BETA3;
    else if (valence > 10)
      return calcBeta(valence);
    else
      return loop_beta[valence-4];
  };

  double calcBeta(int valence) {
    double dval = (double)valence;
    double d = LOOP_MASK_38 + std::cos(2.0 * M_PI / dval) / 4.0;
    double beta = (LOOP_MASK_58 - d * d) / dval;
    return beta;
  };

 private:
  // original mesh
  MeshL* mesh_;

  // subdivided mesh
  MeshL* submesh_;

  // vertices pointer
  std::vector<VertexL*> even_;  // even vertex
  std::vector<VertexL*> odd_;   // odd vertex
};

#endif  // _LOOPSUB_HXX
