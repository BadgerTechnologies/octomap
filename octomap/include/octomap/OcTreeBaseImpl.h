/*
 * OctoMap - An Efficient Probabilistic 3D Mapping Framework Based on Octrees
 * https://octomap.github.io/
 *
 * Copyright (c) 2009-2013, K.M. Wurm and A. Hornung, University of Freiburg
 * All rights reserved.
 * License: New BSD
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Freiburg nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef OCTOMAP_OCTREE_BASE_IMPL_H
#define OCTOMAP_OCTREE_BASE_IMPL_H


#include <list>
#include <limits>
#include <iterator>
#include <stack>
#include <bitset>
#include <functional>

#include <boost/pool/pool.hpp>

#include "octomap_types.h"
#include "OcTreeKey.h"
#include "ScanGraph.h"
#include "OcTreeSpace.h"


namespace octomap {
  
  // forward declaration for NODE children array
  class AbstractOcTreeNode;


  /**
   * OcTree base class, to be used with with any kind of OcTreeDataNode.
   *
   * This tree implementation has a maximum depth limited by key_type
   *
   * \note You should probably not use this class directly, but
   * OcTreeBase or OccupancyOcTreeBase instead
   *
   * \tparam NODE Node class to be used in tree (usually derived from
   *    OcTreeDataNode)
   * \tparam INTERFACE Interface to be derived from, should be either
   *    AbstractOcTree or AbstractOccupancyOcTree
   */
  template <class NODE,class INTERFACE>
  class OcTreeBaseImpl : public INTERFACE, public OcTreeSpace {

  public:
    /// Make the templated NODE type available from the outside
    typedef NODE NodeType;

    // the actual iterator implementation is included here
    // as a member from this file
    #include <octomap/OcTreeIterator.hxx>
    
    OcTreeBaseImpl(double resolution);
    virtual ~OcTreeBaseImpl();

    /// Deep copy constructor
    OcTreeBaseImpl(const OcTreeBaseImpl<NODE,INTERFACE>& rhs);


    /**
     * Swap contents of two octrees, i.e., only the underlying
     * pointer / tree structure. You have to ensure yourself that the
     * metadata (resolution etc) matches. No memory is cleared
     * in this function
     */
    void swapContent(OcTreeBaseImpl<NODE,INTERFACE>& rhs);

    /// Comparison between two octrees, all meta data, all
    /// nodes, and the structure must be identical
    bool operator== (const OcTreeBaseImpl<NODE,INTERFACE>& rhs) const;

    std::string getTreeType() const {return "OcTreeBaseImpl";}

    /// Change the resolution of the octree, scaling all voxels.
    /// This will not preserve the (metric) scale!
    virtual void setResolution(double r);
    virtual double getResolution() const {return OcTreeSpace::getResolution();}

    virtual unsigned int getTreeDepth() const {return OcTreeSpace::getTreeDepth();}
    /// Alter the tree depth. If the tree is not empty, the tree will be
    /// altered to match the new depth, potentially discarding data.
    virtual void setTreeDepth (unsigned int depth);

    /**
     * Clear KeyRay vector to minimize unneeded memory. This is only
     * useful for the StaticMemberInitializer classes, don't call it for
     * an octree that is actually used.
     */
    void clearKeyRays(){
      keyrays.clear();
    }
    
    // -- Tree structure operations formerly contained in the nodes ---
   
    /// Creates (allocates) the i-th child of the node. @return ptr to newly create NODE
    NODE* createNodeChild(NODE* node, unsigned int childIdx);
    
    /// Deletes the i-th child of the node
    void deleteNodeChild(NODE* node, unsigned int childIdx);
    
    /// @return ptr to child number childIdx of node
    NODE* getNodeChild(NODE* node, unsigned int childIdx) const;
    
    /// @return const ptr to child number childIdx of node
    const NODE* getNodeChild(const NODE* node, unsigned int childIdx) const;

    /// Set child number childIdx of node
    void setNodeChild(NODE* node, unsigned int childIdx, NODE* child);

    /// A node is collapsible if all children exist, don't have children of their own
    /// and have the same occupancy value
    virtual bool isNodeCollapsible(const NODE* node) const;
    
    /** 
     * Safe test if node has a child at index childIdx.
     * First tests if there are any children. Replaces node->childExists(...)
     * \return true if the child at childIdx exists
     */
    bool nodeChildExists(const NODE* node, unsigned int childIdx) const;
    
    /** 
     * Safe test if node has any children. Replaces node->hasChildren(...)
     * \return true if node has at least one child
     */
    bool nodeHasChildren(const NODE* node) const;
    
    /**
     * Expands a node (reverse of pruning): All children are created and
     * their occupancy probability is set to the node's value.
     *
     * You need to verify that this is indeed a pruned node (i.e. not a
     * leaf at the lowest level)
     *
     */
    virtual void expandNode(NODE* node);
    
    /**
     * Prunes a node when it is collapsible
     * @return true if pruning was successful
     */
    virtual bool pruneNode(NODE* node);
    
    
    // --------

    /**
     * \return Pointer to the root node of the tree. This pointer
     * should not be modified or deleted externally, the OcTree
     * manages its memory itself. In an empty tree, root is NULL.
     */
    inline NODE* getRoot() const { return root; }

    /** 
     *  Search node at specified depth given a 3d point (depth=0: search full tree depth).
     *  You need to check if the returned node is NULL, since it can be in unknown space.
     *  @return pointer to node if found, NULL otherwise
     */
    NODE* search(double x, double y, double z, unsigned int depth = 0) const;

    /**
     *  Search node at specified depth given a 3d point (depth=0: search full tree depth)
     *  You need to check if the returned node is NULL, since it can be in unknown space.
     *  @return pointer to node if found, NULL otherwise
     */
    NODE* search(const point3d& value, unsigned int depth = 0) const;

    /**
     *  Search a node at specified depth given an addressing key (depth=0: search full tree depth)
     *  You need to check if the returned node is NULL, since it can be in unknown space.
     *  @return pointer to node if found, NULL otherwise
     */
    NODE* search(const OcTreeKey& key, unsigned int depth = 0) const;

    /**
     *  Delete a node (if exists) given a 3d point. Will always
     *  delete at the lowest level unless depth !=0, and expand pruned inner nodes as needed.
     *  Pruned nodes at level "depth" will directly be deleted as a whole.
     */
    bool deleteNode(double x, double y, double z, unsigned int depth = 0);

    /** 
     *  Delete a node (if exists) given a 3d point. Will always
     *  delete at the lowest level unless depth !=0, and expand pruned inner nodes as needed.
     *  Pruned nodes at level "depth" will directly be deleted as a whole.
     */
    bool deleteNode(const point3d& value, unsigned int depth = 0);

    /** 
     *  Delete a node (if exists) given an addressing key. Will always
     *  delete at the lowest level unless depth !=0, and expand pruned inner nodes as needed.
     *  Pruned nodes at level "depth" will directly be deleted as a whole.
     */
    bool deleteNode(const OcTreeKey& key, unsigned int depth = 0);

    using DeletionCallback = std::function<void(const OcTreeBaseImpl<NODE, INTERFACE>* tree, const NODE* node, const OcTreeKey& key, unsigned int depth)>;

    /**
     * Delete all nodes in (or out) of the given AABB.
     *
     * Delete all the nodes in the given axis-aligned bounding box (AABB).
     * This method efficiently descends the tree only visiting the nodes in
     * (or out) of the box. If invert is set, nodes outside the box are
     * deleted (exclusive), otherwise nodes inside or on the box are deleted
     * (inclusive).
     *
     * @param min    Bottom-left corner of box (inclusive).
     * @param max    Top-right corner of box (inclusive).
     * @param invert If true, delete everything out of the box, otherwise delete in the box.
     * @param deletion_notifier If set, a callback to call when a leaf node is
     *               about to be deleted.
     */
    void deleteAABB(const point3d& min, const point3d& max, bool invert = false,
                    DeletionCallback deletion_notifier = DeletionCallback());

    /**
     * Delete all nodes in (or out) of the given AABB.
     *
     * Delete all the nodes in the given axis-aligned bounding box (AABB).
     * This method efficiently descends the tree only visiting the nodes in
     * (or out) of the box. If invert is set, nodes outside the box are
     * deleted (exclusive), otherwise nodes inside or on the box are deleted
     * (inclusive).
     *
     * @param min    Bottom-left corner of box (inclusive).
     * @param max    Top-right corner of box (inclusive).
     * @param invert If true, delete everything out of the box, otherwise delete in the box.
     * @param deletion_notifier If set, a callback to call when a leaf node is
     *               about to be deleted.
     */
    void deleteAABB(const OcTreeKey& min, const OcTreeKey& max, bool invert = false,
                    DeletionCallback deletion_notifier = DeletionCallback());

    /// Deletes the complete tree structure
    void clear();

    /**
     * Lossless compression of the octree: A node will replace all of its eight
     * children if they have identical values. You usually don't have to call
     * prune() after a regular occupancy update, updateNode() incrementally
     * prunes all affected nodes.
     */
    virtual void prune();

    /// Expands all pruned nodes (reverse of prune())
    /// \note This is an expensive operation, especially when the tree is nearly empty!
    virtual void expand();

    // -- statistics  ----------------------

    /// \return The number of nodes in the tree
    virtual inline size_t size() const { return tree_size; }

    /// \return Memory usage of the complete octree in bytes (may vary between architectures)
    virtual size_t memoryUsage() const;

    /// \return Memory usage of a single octree node
    virtual inline size_t memoryUsageNode() const {return sizeof(NODE); }

    /// \return Memory usage of a full grid of the same size as the OcTree in bytes (for comparison)
    /// \note this can be larger than the adressable memory - size_t may not be enough to hold it!
    unsigned long long memoryFullGrid() const;

    double volume();

    /// Size of OcTree (all known space) in meters for x, y and z dimension
    virtual void getMetricSize(double& x, double& y, double& z);
    /// Size of OcTree (all known space) in meters for x, y and z dimension
    virtual void getMetricSize(double& x, double& y, double& z) const;
    /// minimum value of the bounding box of all known space in x, y, z
    virtual void getMetricMin(double& x, double& y, double& z);
    /// minimum value of the bounding box of all known space in x, y, z
    void getMetricMin(double& x, double& y, double& z) const;
    /// maximum value of the bounding box of all known space in x, y, z
    virtual void getMetricMax(double& x, double& y, double& z);
    /// maximum value of the bounding box of all known space in x, y, z
    void getMetricMax(double& x, double& y, double& z) const;

    /// Traverses the tree to calculate the total number of nodes
    size_t calcNumNodes() const;

    /// Traverses the tree to calculate the total number of leaf nodes
    size_t getNumLeafNodes() const;


    // -- access tree nodes  ------------------

    /// return centers of leafs that do NOT exist (but could) in a given bounding box
    void getUnknownLeafCenters(point3d_list& node_centers, point3d pmin, point3d pmax, unsigned int depth = 0) const;


    // -- raytracing  -----------------------

   /**
    * Traces a ray from origin to end (excluding), returning an
    * OcTreeKey of all nodes traversed by the beam. You still need to check
    * if a node at that coordinate exists (e.g. with search()).
    *
    * @param origin start coordinate of ray
    * @param end end coordinate of ray
    * @param ray KeyRay structure that holds the keys of all nodes traversed by the ray, excluding "end"
    * @return Success of operation. Returning false usually means that one of the coordinates is out of the OcTree's range
    */
    bool computeRayKeys(const point3d& origin, const point3d& end, KeyRay& ray) const;


   /**
    * Traces a ray from origin to end (excluding), returning the
    * coordinates of all nodes traversed by the beam. You still need to check
    * if a node at that coordinate exists (e.g. with search()).
    * @note: use the faster computeRayKeys method if possible.
    * 
    * @param origin start coordinate of ray
    * @param end end coordinate of ray
    * @param ray KeyRay structure that holds the keys of all nodes traversed by the ray, excluding "end"
    * @return Success of operation. Returning false usually means that one of the coordinates is out of the OcTree's range
    */
    bool computeRay(const point3d& origin, const point3d& end, std::vector<point3d>& ray);


    // file IO

    /**
     * Read all nodes from the input stream (without file header),
     * for this the tree needs to be already created.
     * For general file IO, you
     * should probably use AbstractOcTree::read() instead.
     */
    std::istream& readData(std::istream &s);

    /// Write complete state of tree to stream (without file header) unmodified.
    /// Pruning the tree first produces smaller files (lossless compression)
    std::ostream& writeData(std::ostream &s) const;

    typedef leaf_iterator iterator;

    /// @return beginning of the tree as leaf iterator
    iterator begin(unsigned char maxDepth=0) const {return iterator(this, maxDepth);}
    /// @return end of the tree as leaf iterator
    const iterator end() const {return leaf_iterator_end;} // TODO: RVE?

    /// @return beginning of the tree as leaf iterator
    leaf_iterator begin_leafs(unsigned char maxDepth=0) const {return leaf_iterator(this, maxDepth);}
    /// @return end of the tree as leaf iterator
    const leaf_iterator end_leafs() const {return leaf_iterator_end;}

    /// @return beginning of the tree as leaf iterator in a bounding box
    leaf_bbx_iterator begin_leafs_bbx(const OcTreeKey& min, const OcTreeKey& max, unsigned char maxDepth=0) const {
      return leaf_bbx_iterator(this, min, max, maxDepth);
    }
    /// @return beginning of the tree as leaf iterator in a bounding box
    leaf_bbx_iterator begin_leafs_bbx(const point3d& min, const point3d& max, unsigned char maxDepth=0) const {
      return leaf_bbx_iterator(this, min, max, maxDepth);
    }
    /// @return end of the tree as leaf iterator in a bounding box
    const leaf_bbx_iterator end_leafs_bbx() const {return leaf_iterator_bbx_end;}

    /// @return beginning of the tree as iterator to all nodes (incl. inner)
    tree_iterator begin_tree(unsigned char maxDepth=0) const {return tree_iterator(this, maxDepth);}
    /// @return end of the tree as iterator to all nodes (incl. inner)
    const tree_iterator end_tree() const {return tree_iterator_end;}

 protected:
    /// Constructor to enable derived classes to change tree constants.
    /// This usually requires a re-implementation of some core tree-traversal functions as well!
    OcTreeBaseImpl(double resolution, unsigned int tree_depth, unsigned int tree_max_val);

    /// initialize non-trivial members, helper for constructors
    void init();

    /// recalculates min and max in x, y, z. Does nothing when tree size didn't change.
    void calcMinMax();

    void calcNumNodesRecurs(NODE* node, size_t& num_nodes) const;
    
    /// recursive call of readData()
    std::istream& readNodesRecurs(NODE*, std::istream &s);
    
    /// recursive call of writeData()
    std::ostream& writeNodesRecurs(const NODE*, std::ostream &s) const;

    /// Recusively clone a node and all children.
    NODE* cloneNodeRecurs(NODE* node);

    /// Recursively delete a node and all children. Deallocates memory
    /// but does NOT set the node ptr to NULL.
    void deleteNodeRecurs(NODE* node);

    /// recursive call of deleteNode()
    bool deleteNodeRecurs(NODE* node, unsigned int depth, unsigned int max_depth, const OcTreeKey& key);

    /// Recursively delete a node and all children, calling the deletion_notifier for every leaf node deleted.
    void deleteNodeRecurs(NODE* node,
                          const OcTreeKey& key,
                          unsigned int depth,
                          const DeletionCallback& deletion_notifier);

    /// recursive call of deleteAABB(). Returns true if node has been deleted.
    bool deleteAABBRecurs(const OcTreeKey& min, const OcTreeKey& max,
                          NODE* node, const OcTreeKey& key, unsigned int depth,
                          unsigned int max_depth, bool invert,
                          const DeletionCallback& deletion_notifier);

    /// recursive call of prune()
    void pruneRecurs(NODE* node, unsigned int depth, unsigned int max_depth, unsigned int& num_pruned);

    /// recursive call of expand()
    void expandRecurs(NODE* node, unsigned int depth, unsigned int max_depth);
    
    size_t getNumLeafNodesRecurs(const NODE* parent) const;

  private:
    /// Assignment operator is private: don't (re-)assign octrees
    /// (const-parameters can't be changed) -  use the copy constructor instead.
    OcTreeBaseImpl<NODE,INTERFACE>& operator=(const OcTreeBaseImpl<NODE,INTERFACE>&);

  protected:
    NODE* allocNode();
    void freeNode(NODE* node);
    void allocNodeChildren(NODE* node);
    void deleteNodeChildren(NODE* node);
    void deleteNodeChildren(NODE* node,
                            const OcTreeKey& key,
                            unsigned int depth,
                            const DeletionCallback& deletion_notifier);
    void deleteNodeChildrenIfNecessary(NODE* node);

    NODE* root; ///< Pointer to the root NODE, NULL for empty tree

    size_t tree_size; ///< number of nodes in tree
    /// flag to denote whether the octree extent changed (for lazy min/max eval)
    bool size_changed;

    double max_value[3]; ///< max in x, y, z
    double min_value[3]; ///< min in x, y, z

    /// data structure for ray casting, array for multithreading
    std::vector<KeyRay> keyrays;

    const leaf_iterator leaf_iterator_end;
    const leaf_bbx_iterator leaf_iterator_bbx_end;
    const tree_iterator tree_iterator_end;

    boost::pool<> node_pool;
    boost::pool<> children_pool;
  };

}

#include <octomap/OcTreeBaseImpl.hxx>

#endif
