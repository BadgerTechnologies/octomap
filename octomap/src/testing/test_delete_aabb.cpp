#include <memory>
#include <functional>
#include <octomap/octomap.h>
#include "octomap/AbstractOccupancyOcTree.h"
#include "testing.h"

using namespace std;
using namespace octomap;
using namespace octomath;

int main(int argc, char** argv) {

  float res = 0.01f;
  OcTree tree(res);
  OcTree deleted_tree(res);
  OcTree deleted_tree2(res);
  shared_ptr<OcTree> treep, saved_tree;

  std::cout << "\nDeleting AABBs \n===============================\n";

  // Start by verifying coordToKeyClamped
  const key_type min_key_val = std::numeric_limits<key_type>::min();
  const key_type max_key_val = std::numeric_limits<key_type>::max();
  const key_type avg_key_val = (max_key_val - min_key_val) / 2;
  EXPECT_EQ(avg_key_val, tree.coordToKey(tree.keyToCoord(avg_key_val)));
  EXPECT_EQ(min_key_val, tree.coordToKey(tree.keyToCoord(min_key_val)));
  EXPECT_EQ(max_key_val, tree.coordToKey(tree.keyToCoord(max_key_val)));
  OcTreeKey key;
  tree.coordToKeyClamped(-1e37, -1e37, -1e37, key);
  EXPECT_EQ(min_key_val, key[0]);
  EXPECT_EQ(min_key_val, key[1]);
  EXPECT_EQ(min_key_val, key[2]);
  tree.coordToKeyClamped(1e37, 1e37, 1e37, key);
  EXPECT_EQ(max_key_val, key[0]);
  EXPECT_EQ(max_key_val, key[1]);
  EXPECT_EQ(max_key_val, key[2]);
  tree.coordToKeyClamped(-1e-37, -1e-37, -1e-37, key);
  EXPECT_EQ(avg_key_val, key[0]);
  EXPECT_EQ(avg_key_val, key[1]);
  EXPECT_EQ(avg_key_val, key[2]);

  for(int i=0; i<16; ++i) {
    for(int j=0; j<16; ++j) {
      for(int k=0; k<16; ++k) {
        float value = -1.0;
        if (i >= 1 && i <= 3 && j >= 1 && j <= 3 && k >= 1 && k <= 3) {
          value = 1.0;
        }
        OcTreeKey key(32768+i, 32768+j, 32768+k);
        tree.setNodeValue(key, value);
      }
    }
  }

  EXPECT_EQ(tree.size(), 93);
  treep.reset(new OcTree(tree));
  EXPECT_EQ(treep->size(), 93);
  treep->deleteAABB(point3d(0.001, 0.001, 0.001), point3d(0.039, 0.039, 0.039), false,
                    std::bind([](OcTree* tree, const OcTreeNode* node, const OcTreeKey& key, unsigned int depth)
                              {tree->setNodeValueAtDepth(key, depth, node->getValue());},
                              &deleted_tree,
                              std::placeholders::_2,
                              std::placeholders::_3,
                              std::placeholders::_4));
  EXPECT_EQ(treep->size(), 28);
  saved_tree = treep;
  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(0.001, 0.001, 0.001), point3d(0.039, 0.039, 0.039), true,
                    std::bind([](OcTree* tree, const OcTreeNode* node, const OcTreeKey& key, unsigned int depth)
                              {tree->setNodeValueAtDepth(key, depth, node->getValue());},
                              &deleted_tree2,
                              std::placeholders::_2,
                              std::placeholders::_3,
                              std::placeholders::_4));
  EXPECT_EQ(treep->size(), 79);
  // What is deleted here should exactly match what was not deleted in the
  // same call above when invert was false.
  EXPECT_TRUE(*saved_tree == deleted_tree2);
  EXPECT_TRUE(*treep == deleted_tree);

  tree.clear();
  tree.setNodeValueAtDepth(OcTreeKey(32768+16, 32768+16, 32768+16), 11, 1.0);
  EXPECT_EQ(tree.size(), 12);

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(0.001, 0.001, 0.001), point3d(0.321, 0.321, 0.321));
  EXPECT_EQ(treep->size(), 0);

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(-1e37, -1e37, -1e37), point3d(1e37, 1e37, 1e37));
  EXPECT_EQ(treep->size(), 0);

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(0.001, 0.001, 0.001), point3d(0.319, 0.319, 0.319));
  EXPECT_EQ(treep->size(), 0);

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(0.011, 0.011, 0.011), point3d(0.319, 0.319, 0.319));
  EXPECT_EQ(treep->size(), 3923);

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(0.011, 0.011, 0.011), point3d(30.319, 30.319, 0.319));
  EXPECT_EQ(treep->size(), 3923);

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(-300.011, -30.011, -0.011), point3d(0.309, 0.309, 0.309));
  EXPECT_EQ(treep->size(), 3923);

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(0.011, 0.011, 0.011), point3d(0.309, 0.309, 0.309));
  EXPECT_EQ(treep->size(), 7492);
  saved_tree = treep;

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(0.011, 0.011, 0.011), point3d(0.319, 0.319, 0.319), true);
  EXPECT_EQ(treep->size(), 4515);

  treep.reset(new OcTree(tree));
  deleted_tree.clear();
  treep->deleteAABB(point3d(0.011, 0.011, 0.011), point3d(0.309, 0.309, 0.309), true,
                    std::bind([](OcTree* tree, const OcTreeNode* node, const OcTreeKey& key, unsigned int depth)
                              {tree->setNodeValueAtDepth(key, depth, node->getValue());},
                              &deleted_tree,
                              std::placeholders::_2,
                              std::placeholders::_3,
                              std::placeholders::_4));
  EXPECT_EQ(treep->size(), 7948);
  // What is deleted here should exactly match what was not deleted in the
  // same call above when invert was false.
  EXPECT_TRUE(*saved_tree == deleted_tree);

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(0.071, 0.071, 0.071), point3d(0.079, 0.079, 0.079));
  EXPECT_EQ(treep->size(), 51);

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(0.071, 0.071, 0.071), point3d(0.079, 0.079, 0.079), true);
  EXPECT_EQ(treep->size(), 17);

  treep.reset(new OcTree(tree));
  treep->deleteAABB(point3d(0.071, 0.071, 0.071), point3d(0.059, 0.079, 0.079));
  EXPECT_EQ(treep->size(), 12);

  treep.reset(new OcTree(tree));
  deleted_tree.clear();
  treep->deleteAABB(point3d(0.071, 0.071, 0.071), point3d(0.059, 0.079, 0.079), true,
                    std::bind([](OcTree* tree, const OcTreeNode* node, const OcTreeKey& key, unsigned int depth)
                              {tree->setNodeValueAtDepth(key, depth, node->getValue());},
                              &deleted_tree,
                              std::placeholders::_2,
                              std::placeholders::_3,
                              std::placeholders::_4));
  EXPECT_EQ(treep->size(), 0);
  EXPECT_TRUE(deleted_tree == tree);

  std::cerr << "Test successful.\n";
  return 0;
}
