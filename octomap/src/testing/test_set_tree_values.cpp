#include <memory>
#include <octomap/octomap.h>
#include "testing.h"

using namespace std;
using namespace octomap;
using namespace octomath;

int main(int argc, char** argv) {
    float res = 0.01f;
    OcTree value_tree(res);
    OcTree value_tree2(res);
    OcTree bounds_tree(res);
    OcTree tree(res);
    OcTree expected_tree(res);
    shared_ptr<OcTree> treep;
    OcTreeNode* node;

    std::cout << "\nSetting Tree Values from Other Trees\n===============================\n";

    // First, test using empty trees
    tree.setTreeValues(NULL, false, false);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(NULL, true, false);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(NULL, true, true);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, true);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, true, true);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, &bounds_tree);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, &bounds_tree, true);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, &bounds_tree, true, true);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    // Add a node to the value tree but use an empty bounds tree
    value_tree.setNodeValue(OcTreeKey(), 1.0);
    tree.setTreeValues(&value_tree, &bounds_tree);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, &bounds_tree, true);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, &bounds_tree, true, true);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());

    // Test the case of an empty tree being set to the universe (pruned node at top).
    tree.clear();
    value_tree.clear();
    value_tree.setNodeValueAtDepth(OcTreeKey(), 0, value_tree.getClampingThresMaxLog());
    EXPECT_EQ(tree.size(), 0);
    tree.setTreeValues(&value_tree);
    EXPECT_EQ(tree.size(), 1);
    tree.clear();
    value_tree.clear();

    // Now, test with one leaf in our tree
    point3d singlePt(-0.05f, -0.02f, 1.0f);
    OcTreeKey singleKey, nextKey;
    tree.coordToKeyChecked(singlePt, singleKey);
    tree.updateNode(singleKey, true);
    tree.setTreeValues(&value_tree);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, true);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, false, true);
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.updateNode(singleKey, true);
    // Since the bounds tree is empty, no change should happen
    tree.setTreeValues(&value_tree, &bounds_tree, true, true);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    bounds_tree.updateNode(singleKey, true);
    tree.setTreeValues(&value_tree, &bounds_tree, false, false);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, &bounds_tree, true, false);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, &bounds_tree, false, true);
    // Bounds tree has our key, our tree should be empty now
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());

    // Set the bounds tree to everything.
    bounds_tree.setNodeValueAtDepth(OcTreeKey(), 0, bounds_tree.getClampingThresMax());
    EXPECT_EQ(bounds_tree.size(), 1);
    tree.clear();
    tree.updateNode(singleKey, true);
    tree.setTreeValues(&value_tree, &bounds_tree, false, false);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, &bounds_tree, true, false);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    tree.setTreeValues(&value_tree, &bounds_tree, false, true);
    // Bounds tree has our key, our tree should be empty now
    EXPECT_EQ(tree.size(), 0);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());

    // Now put a single node in the value tree
    value_tree.setNodeValue(singleKey, -1.0);
    tree.setNodeValue(singleKey, 1.0);
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(1.0, node->getLogOdds());
    tree.setTreeValues(&value_tree, true, false);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(1.0, node->getLogOdds());
    tree.setTreeValues(&value_tree, false, false);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(-1.0, node->getLogOdds());
    tree.setNodeValue(singleKey, 1.0);
    tree.setTreeValues(&value_tree, false, true);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(-1.0, node->getLogOdds());

    // Now try the same node in the bounds tree.
    bounds_tree.clear();
    bounds_tree.setNodeValue(singleKey, 1.0);
    tree.setNodeValue(singleKey, 1.0);
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(1.0, node->getLogOdds());
    tree.setTreeValues(&value_tree, &bounds_tree, true, false);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(1.0, node->getLogOdds());
    tree.setTreeValues(&value_tree, &bounds_tree, false, false);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(-1.0, node->getLogOdds());
    tree.setNodeValue(singleKey, 1.0);
    tree.setTreeValues(&value_tree, &bounds_tree, false, true);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(-1.0, node->getLogOdds());

    // Set the bounds tree to everything.
    bounds_tree.setNodeValueAtDepth(OcTreeKey(), 0, bounds_tree.getClampingThresMax());
    EXPECT_EQ(bounds_tree.size(), 1);
    tree.setNodeValue(singleKey, 1.0);
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(1.0, node->getLogOdds());
    tree.setTreeValues(&value_tree, &bounds_tree, true, false);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(1.0, node->getLogOdds());
    tree.setTreeValues(&value_tree, &bounds_tree, false, false);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(-1.0, node->getLogOdds());
    tree.setNodeValue(singleKey, 1.0);
    tree.setTreeValues(&value_tree, &bounds_tree, false, true);
    EXPECT_EQ(tree.size(), 17);
    EXPECT_EQ(tree.size(), tree.calcNumNodes());
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(-1.0, node->getLogOdds());

    // Test having a pruned inner node for a value tree when using maximum.
    tree.setNodeValue(singleKey, 1.0);
    value_tree.setNodeValueAtDepth(singleKey, 15, -1.0);
    EXPECT_EQ(value_tree.size(), 16);
    EXPECT_EQ(tree.size(), 17);
    tree.setTreeValues(&value_tree, &bounds_tree, true, false);
    EXPECT_EQ(tree.size(), 24);
    node = tree.search(singleKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(1.0, node->getLogOdds());
    nextKey = singleKey;
    nextKey[1] += 1;
    node = tree.search(nextKey);
    EXPECT_TRUE(node);
    EXPECT_EQ(-1.0, node->getLogOdds());

    // Test no value tree with a non-overlapping bounds tree with delete first set.
    tree.clear();
    bounds_tree.clear();
    nextKey[1] += 2;
    tree.setNodeValue(singleKey, 1.0);
    bounds_tree.setNodeValue(nextKey, 1.0);
    EXPECT_EQ(17, tree.size());
    EXPECT_EQ(17, bounds_tree.size());
    treep.reset(new OcTree(tree));
    treep->setTreeValues(NULL, &bounds_tree, false, true);
    EXPECT_EQ(17, tree.size());
    EXPECT_TRUE(tree == *treep);

    // Test delete first with an empty value tree, a complex bounds tree and
    // value tree with the value tree completely inside (result should be
    // empty).
    tree.clear();
    bounds_tree.clear();
    for(int i=-10; i<=10; ++i) {
      for(int j=-10; j<=10; ++j) {
        for(int k=-10; k<=10; ++k) {
          OcTreeKey key(32768+i, 32768+j, 32768+k);
          bounds_tree.setNodeValue(key, 1.0);
        }
      }
    }
    for(int i=-7; i<=9; ++i) {
      for(int j=-5; j<=3; ++j) {
        for(int k=-2; k<=8; ++k) {
          OcTreeKey key(32768+i, 32768+j, 32768+k);
          tree.setNodeValue(key, 1.0);
        }
      }
    }
    tree.setTreeValues(NULL, &bounds_tree, false, true);
    EXPECT_EQ(0, tree.size());

    // Now, make more complex merging scenarios.
    tree.clear();
    value_tree.clear();
    bounds_tree.setNodeValueAtDepth(OcTreeKey(), 0, bounds_tree.getClampingThresMax());
    for(int i=0; i<4; ++i) {
      for(int j=0; j<4; ++j) {
        for(int k=0; k<4; ++k) {
          float value1 = -1.0;
          float value2 = 1.0;
          if (i >= 1 && i <= 2 && j >= 1 && j <= 2 && k >= 1 && k <= 2) {
            value1 = 1.0;
            value2 = -1.0;
          }
          OcTreeKey key(32768+i, 32768+j, 32768+k);
          tree.setNodeValue(key, value1);
          value_tree.setNodeValue(key, value2);
        }
      }
    }
    OcTreeKey search_key(32768+2, 32768+2, 32768+2);
    EXPECT_EQ(4*4*4 + 4*4*4/8 + 1*15, tree.size());
    EXPECT_EQ(4*4*4 + 4*4*4/8 + 1*15, value_tree.size());
    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, &bounds_tree, true, false);
    EXPECT_EQ(15, treep->size());
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    node = treep->search(search_key);
    EXPECT_EQ(1.0, node->getLogOdds());
    expected_tree.clear();
    expected_tree.setNodeValueAtDepth(OcTreeKey(32768+2, 32768+2, 32768+2), 14, 1.0);
    EXPECT_TRUE(expected_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, &bounds_tree, false, false);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(value_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, &bounds_tree, false, true);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(value_tree == *treep);

    // Now try with bounds limited.
    expected_tree.clear();
    treep.reset(new OcTree(tree));
    expected_tree.swapContent(*treep);
    expected_tree.setNodeValueAtDepth(OcTreeKey(32768+1, 32768+1, 32768+1), 15, 1.0);
    bounds_tree.clear();
    bounds_tree.setNodeValueAtDepth(OcTreeKey(32768+1, 32768+1, 32768+1), 15, 1.0);
    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, &bounds_tree, true, false);
    EXPECT_EQ(4*4*4 - 8 + 4*4*4/8 + 1*15, treep->size());
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(expected_tree == *treep);

    expected_tree.setNodeValue(OcTreeKey(32768+1, 32768+1, 32768+1), -1.0);
    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, &bounds_tree, false, false);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(expected_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, &bounds_tree, false, true);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(expected_tree == *treep);

    tree.clear();
    value_tree.clear();
    value_tree2.clear();
    bounds_tree.clear();

    for(int i=-4; i<4; ++i) {
      for(int j=-4; j<4; ++j) {
        for(int k=-4; k<4; ++k) {
          float value1 = -200.0;
          float value2 = -200.0;
          float value3 = -200.0;
          if (i >= -4 && i <= -3 && j >= -4 && j <= -3 && k >= -4 && k <= -3)
          {
            value1 = -1.0;
          }
          else if( i >= -4 && i < 0 && j >= -4 && j < 0 && k >= -4 && k < 0)
          {
            value1 = 1.0;
          }
          else if(i == 0 && j == 0 && k == 0)
          {
            value1 = 1.0;
          }
          else if(i == 2 && j == 2 && k == 2)
          {
            value1 = 1.0;
          }
          else if (i >= 2 && i <= 3 && j >= 2 && j <= 3 && k >= 2 && k <= 3)
          {
            value1 = -1.0;
          }
          if( i >= -4 && i < 0 && j >= -4 && j < 0 && k >= -4 && k < 0)
          {
            value2 = 1.0;
          }
          else if( i >= 1 && i < 4 && j >= 1 && j < 4 && k >= 1 && k < 4)
          {
            value2 = -1.0;
          }
          if(i == -2 && j == -2 && k == -2)
          {
            value3 = -1.0;
          }
          else if(i == 0 && j == 0 && k == 0)
          {
            value3 = -1.0;
          }
          else if(i == 1 && j == 1 && k == 1)
          {
            value3 = 1.0;
          }
          else if (i >= 2 && i <= 3 && j >= 2 && j <= 3 && k >= 2 && k <= 3)
          {
            value3 = 1.0;
          }

          OcTreeKey key(32768+i, 32768+j, 32768+k);
          if (value1 > -100.0)
          {
            tree.setNodeValue(key, value1);
          }
          if (value2 > -100.0)
          {
            value_tree.setNodeValue(key, value2);
          }
          if (value3 > -100.0)
          {
            value_tree2.setNodeValue(key, value3);
          }
        }
      }
    }
    bounds_tree.setNodeValueAtDepth(OcTreeKey(32768-2, 32768-2, 32768-2), 14, 1.0);
    bounds_tree.setNodeValueAtDepth(OcTreeKey(32768+3, 32768+3, 32768+3), 15, 1.0);
    treep.reset(new OcTree(tree));
    EXPECT_EQ((1 + 8) + 2 + 15 + 8 + 14, treep->size());
    EXPECT_EQ(14 + 4*3+2*3+1 + 8 + 15, value_tree.size());
    EXPECT_EQ(16 + 2 + 2 + 15, value_tree2.size());

    treep->setTreeValues(&value_tree, &bounds_tree, true, false);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2, &bounds_tree, true, false);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_EQ(15 + 17, treep->size());
    expected_tree.clear();
    expected_tree.setNodeValueAtDepth(OcTreeKey(32768-2, 32768-2, 32768-2), 14, 1.0);
    expected_tree.setNodeValue(OcTreeKey(32768, 32768, 32768), 1.0);
    expected_tree.setNodeValueAtDepth(OcTreeKey(32768+3, 32768+3, 32768+3), 15, 1.0);
    EXPECT_TRUE(expected_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, true, false);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2, true, false);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    for(int i=1; i<4; ++i) {
      for(int j=1; j<4; ++j) {
        for(int k=1; k<4; ++k) {
          float v=-200.0;
          if (i==1 && j==1 && k==1) {
            v=1.0;
          } else if (i==1 || j==1 || k==1) {
            v=-1.0;
          }
          if (v > -100.0) {
            expected_tree.setNodeValue(OcTreeKey(32768+i, 32768+j, 32768+k), v);
          }
        }
      }
    }
    EXPECT_TRUE(expected_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, true);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2, true);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    expected_tree.swapContent(*treep);
    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree2, true);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree, true);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(expected_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    expected_tree.swapContent(*treep);
    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree2);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(expected_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, &bounds_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2, &bounds_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    expected_tree.swapContent(*treep);
    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree2, &bounds_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree, &bounds_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2, &bounds_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(expected_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, &value_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2, &value_tree2);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    expected_tree.swapContent(*treep);
    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree2, &value_tree2);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree, &value_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2, &value_tree2);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(expected_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, &value_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2, &value_tree2);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    expected_tree.swapContent(*treep);
    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree2, &value_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree, &value_tree2);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    treep->setTreeValues(&value_tree2, &value_tree);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_FALSE(expected_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, false, true);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(value_tree == *treep);
    treep->setTreeValues(&value_tree2, false, true);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_TRUE(value_tree2 == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree, &bounds_tree, false, true);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_FALSE(value_tree == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&value_tree2, &bounds_tree, false, true);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());
    EXPECT_FALSE(value_tree2 == *treep);

    treep.reset(new OcTree(tree));
    treep->setTreeValues(NULL, &tree, false, true);
    EXPECT_EQ(treep->size(), 0);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());

    treep.reset(new OcTree(tree));
    treep->setTreeValues(&tree, false, false,
                         [](const OcTree::NodeType*, OcTree::NodeType* node){node->setLogOdds(1.0);});
    EXPECT_EQ(treep->size(), 32);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());

    treep->clear();
    float max_log = tree.getClampingThresMaxLog();
    // Ugly, but in C++11 you can't capture in a lambda and use it as a function pointer.
    // Instead, use bind to send the value to set to the lambda
    treep->setTreeValues(&tree, false, false,
                         std::bind([](OcTree::NodeType* node, float value){node->setLogOdds(value);},
                                   std::placeholders::_2, max_log));
    EXPECT_EQ(treep->size(), 32);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());

    treep->clear();
    treep->setTreeValues(&tree, false, false,
                         [](const OcTree::NodeType*, OcTree::NodeType* node){node->setLogOdds(1.0);});
    treep->setTreeValues(&value_tree, false, false,
                         [](const OcTree::NodeType*, OcTree::NodeType* node){node->setLogOdds(1.0);});
    treep->setTreeValues(&value_tree2, false, false,
                         [](const OcTree::NodeType*, OcTree::NodeType* node){node->setLogOdds(1.0);});
    EXPECT_EQ(treep->size(), 57);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());

    tree.clear();
    value_tree.clear();
    value_tree2.clear();
    for(int i=0; i<16; ++i) {
      for(int j=0; j<16; ++j) {
        for(int k=0; k<16; ++k) {
          float value1 = -200.0;
          float value2 = -200.0;
          float value3 = -200.0;
          if(i<2 || k>13 ) {
            value1 = -1.0;
          }
          else if(i<5 || j>3) {
            value2 = .5;
          }
          else {
            value3 = 1.0;
          }
          OcTreeKey key(32768+i, 32768+j, 32768+k);
          if (value1 > -100.0)
          {
            tree.setNodeValue(key, value1);
          }
          if (value2 > -100.0)
          {
            value_tree.setNodeValue(key, value2);
          }
          if (value3 > -100.0)
          {
            value_tree2.setNodeValue(key, value3);
          }
        }
      }
    }

    treep->clear();
    treep->setTreeValues(&tree, false, false,
                         [](const OcTree::NodeType*, OcTree::NodeType* node){node->setLogOdds(1.0);});
    treep->setTreeValues(&value_tree, false, false,
                         [](const OcTree::NodeType*, OcTree::NodeType* node){node->setLogOdds(1.0);});
    treep->setTreeValues(&value_tree2, false, false,
                         [](const OcTree::NodeType*, OcTree::NodeType* node){node->setLogOdds(1.0);});
    EXPECT_EQ(treep->size(), 13);
    EXPECT_EQ(treep->size(), treep->calcNumNodes());

    std::cerr << "Test successful.\n";
    return 0;

}
