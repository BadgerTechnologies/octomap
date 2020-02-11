#include <stdio.h>
#include <string>
#ifdef _WIN32
  #include <Windows.h>  // to define Sleep()
#else
  #include <unistd.h>   // POSIX sleep()
#endif


#include <octomap/octomap.h>
#include <octomap/OcTreeStamped.h>
#include <octomap/math/Utils.h>
#include "testing.h"
 
using namespace std;
using namespace octomap;
using namespace octomath;

int main(int argc, char** argv) {

  if (argc != 2){
    std::cerr << "Error: you need to specify a test as argument" << std::endl;
    return 1; // exit 1 means failure
  }
  std::string test_name (argv[1]);


  // ------------------------------------------------------------
  if (test_name == "MathVector") {
    // test constructors
    Vector3* twos = new Vector3();        
    Vector3* ones = new Vector3(1,1,1);    
    for (int i=0;i<3;i++) {
      (*twos)(i) = 2;
    }  
    // test basic operations
    Vector3 subtraction = *twos - *ones;
    Vector3 addition = *twos + *ones;
    Vector3 multiplication = *twos * 2.;
  
    for (int i=0;i<3;i++) {
      EXPECT_FLOAT_EQ (subtraction(i), 1.);
      EXPECT_FLOAT_EQ (addition(i), 3.);
      EXPECT_FLOAT_EQ (multiplication(i), 4.);
    }

    // copy constructor
    Vector3 rotation =  *ones;

    // rotation
    rotation.rotate_IP (M_PI, 1., 0.1);
    EXPECT_FLOAT_EQ (rotation.x(), 1.2750367);
    EXPECT_FLOAT_EQ (rotation.y(), (-1.1329513));
    EXPECT_FLOAT_EQ (rotation.z(), 0.30116868);
  
  // ------------------------------------------------------------
  } else if (test_name == "MathPose") {
    // constructors  
    Pose6D a (1.0, 0.1, 0.1, 0.0, 0.1, M_PI/4. );
    Pose6D b;

    Vector3 trans(1.0, 0.1, 0.1);
    Quaternion rot(0.0, 0.1, M_PI/4.);
    Pose6D c(trans, rot);

    // comparator
    EXPECT_TRUE ( a == c);
    // toEuler
    EXPECT_FLOAT_EQ (c.yaw() , M_PI/4.);

    // transform
    Vector3 t = c.transform (trans);
    EXPECT_FLOAT_EQ (t.x() , 1.6399229);
    EXPECT_FLOAT_EQ (t.y() , 0.8813442);
    EXPECT_FLOAT_EQ (t.z() , 0.099667005);

    // inverse transform
    Pose6D c_inv = c.inv();
    Vector3 t2 = c_inv.transform (t);
    EXPECT_FLOAT_EQ (t2.x() , trans.x());
    EXPECT_FLOAT_EQ (t2.y() , trans.y());
    EXPECT_FLOAT_EQ (t2.z() , trans.z());

  // ------------------------------------------------------------
  } else if (test_name == "InsertRay") {
    double p = 0.5;
    EXPECT_FLOAT_EQ(p, probability(logodds(p)));
    p = 0.1;
    EXPECT_FLOAT_EQ(p, probability(logodds(p)));
    p = 0.99;
    EXPECT_FLOAT_EQ(p, probability(logodds(p)));

    float l = 0;
    EXPECT_FLOAT_EQ(l, logodds(probability(l)));
    l = -4;
    EXPECT_FLOAT_EQ(l, logodds(probability(l)));
    l = 2;
    EXPECT_FLOAT_EQ(l, logodds(probability(l)));


    OcTree tree (0.05);
    tree.setProbHit(0.7);
    tree.setProbMiss(0.4);

    point3d origin (0.01, 0.01, 0.02);
    point3d point_on_surface (2.01,0.01,0.01);
  
    for (int i=0; i<360; i++) {    
      for (int j=0; j<360; j++) {
        EXPECT_TRUE (tree.insertRay(origin, origin+point_on_surface));
        point_on_surface.rotate_IP (0,0,DEG2RAD(1.));
      }
      point_on_surface.rotate_IP (0,DEG2RAD(1.),0);
    }
    EXPECT_TRUE (tree.writeBinary("sphere_rays.bt"));
    EXPECT_EQ (tree.size(), 50493 + 8 * tree.getTreeDepth());
  
  // ------------------------------------------------------------
  // ray casting is now in "test_raycasting.cpp"

  // ------------------------------------------------------------
  // insert scan test
  // insert graph node test
  // write graph test
  } else if (test_name == "InsertScan") {
    Pointcloud* measurement = new Pointcloud();
  
    point3d origin (0.01, 0.01, 0.02);
    point3d point_on_surface (2.01, 0.01, 0.01);
  
    for (int i=0; i<360; i++) {
      for (int j=0; j<360; j++) {
        point3d p = origin+point_on_surface;
        measurement->push_back(p);
        point_on_surface.rotate_IP (0,0,DEG2RAD(1.));
      }
      point_on_surface.rotate_IP (0,DEG2RAD(1.),0);
    }
  
    OcTree tree (0.05);
    tree.insertPointCloud(*measurement, origin);
    EXPECT_EQ (tree.size(), 53837 + 8 * tree.getTreeDepth());

    ScanGraph* graph = new ScanGraph();
    Pose6D node_pose (origin.x(), origin.y(), origin.z(),0.0,0.0,0.0);
    graph->addNode(measurement, node_pose);
    EXPECT_TRUE (graph->writeBinary("test.graph"));
    delete graph;
  // ------------------------------------------------------------
  // graph read file test
  } else if (test_name == "ReadGraph") {
    // not really meaningful, see better test in "test_scans.cpp"
    ScanGraph graph;
    EXPECT_TRUE (graph.readBinary("test.graph"));
  // ------------------------------------------------------------

  } else if (test_name == "StampedTree") {
    OcTreeStamped stamped_tree (0.05);
    // fill tree
    for (int x=-20; x<20; x++) 
      for (int y=-20; y<20; y++) 
        for (int z=-20; z<20; z++) {
          point3d p (x*0.05+0.01, y*0.05+0.01, z*0.05+0.01);
          stamped_tree.updateNode(p, true); // integrate 'occupied' measurement 
        }
    // test if update times set
    point3d query (0.1, 0.1, 0.1);
    OcTreeNodeStamped* result = stamped_tree.search (query);
    EXPECT_TRUE (result);
    unsigned int tree_time = stamped_tree.getLastUpdateTime();
    unsigned int node_time = result->getTimestamp();
    std::cout << "After 1st update (cube): Tree time " <<tree_time << "; node(0.1, 0.1, 0.1) time " << result->getTimestamp() << std::endl;
    EXPECT_TRUE (tree_time > 0);
    #ifdef _WIN32
      Sleep(1000);
    #else
      sleep(1);
    #endif
    stamped_tree.integrateMissNoTime(result);  // reduce occupancy, no time update
    std::cout << "After 2nd update (single miss): Tree time " <<tree_time << "; node(0.1, 0.1, 0.1) time " << node_time << std::endl;
    EXPECT_EQ  (node_time, result->getTimestamp()); // node time updated?
    point3d query2 = point3d  (0.1, 0.1, 0.3);
    stamped_tree.updateNode(query2, true); // integrate 'occupied' measurement
    OcTreeNodeStamped* result2 = stamped_tree.search (query2);
    EXPECT_TRUE (result2);
    result = stamped_tree.search (query);
    EXPECT_TRUE (result);
    std::cout << "After 3rd update (single hit at (0.1, 0.1, 0.3): Tree time " << stamped_tree.getLastUpdateTime() << "; node(0.1, 0.1, 0.1) time " << result->getTimestamp()
        << "; node(0.1, 0.1, 0.3) time " << result2->getTimestamp() << std::endl;
    EXPECT_TRUE (result->getTimestamp() < result2->getTimestamp()); // result2 has been updated
    EXPECT_EQ(result2->getTimestamp(), stamped_tree.getLastUpdateTime());
  // ------------------------------------------------------------
  } else if (test_name == "OcTreeKey") {
    OcTree tree (0.05);  
    point3d p(0.0,0.0,0.0);
    OcTreeKey key;
    tree.coordToKeyChecked(p, key);
    point3d p_inv = tree.keyToCoord(key);
    EXPECT_FLOAT_EQ (0.025, p_inv.x());
    EXPECT_FLOAT_EQ (0.025, p_inv.y());
    EXPECT_FLOAT_EQ (0.025, p_inv.z());

  // ------------------------------------------------------------
  } else if (test_name == "setTreeDepth") {
    OcTree tree (0.05), tree2 (0.05);
    point3d p(0.0,0.0,0.0);
    OcTreeKey key;
    tree.coordToKeyChecked(p, key);
    tree.setNodeValueAtDepth(key, KEY_BIT_WIDTH/2 + 1, 1.0);
    tree2.setNodeValueAtDepth(key, KEY_BIT_WIDTH/2 + 1, 1.0);
    EXPECT_EQ(tree.size(), KEY_BIT_WIDTH/2 + 2);
    EXPECT_EQ(tree.search(key, 1)->getValue(), 1.0);
    EXPECT_EQ(tree.getTreeLogOdds(), 1.0);
    tree.setTreeDepth(KEY_BIT_WIDTH/2);
    EXPECT_EQ(tree.size(), 2);
    tree.coordToKeyChecked(p, key);
    EXPECT_EQ(tree.search(key, 1)->getValue(), 1.0);
    EXPECT_EQ(tree.getTreeLogOdds(), 1.0);
    tree.setTreeDepth(KEY_BIT_WIDTH);
    EXPECT_EQ(tree.size(), KEY_BIT_WIDTH/2 + 2);
    tree.coordToKeyChecked(p, key);
    EXPECT_EQ(tree.search(key, 1)->getValue(), 1.0);
    EXPECT_EQ(tree.getTreeLogOdds(), 1.0);
    EXPECT_TRUE(tree == tree2);
    tree.setNodeValueAtDepth(key, 0, 1.0);
    EXPECT_EQ(tree.size(), 1);
    tree.setTreeDepth(KEY_BIT_WIDTH-1);
    EXPECT_EQ(tree.size(), 1);
    tree.setTreeDepth(KEY_BIT_WIDTH/2);
    EXPECT_EQ(tree.size(), 1);
    tree.setTreeDepth(5);
    EXPECT_EQ(tree.size(), 1);
    tree.setTreeDepth(KEY_BIT_WIDTH);
    EXPECT_EQ(tree.size(), 1 + 8*(KEY_BIT_WIDTH-5+1));
    tree.clear();
    tree.coordToKeyChecked(p, key);
    tree.setNodeValueAtDepth(key, KEY_BIT_WIDTH-3, 3.0);
    key[0] -= 1;
    tree.setNodeValueAtDepth(key, KEY_BIT_WIDTH-2, 2.0);
    key[2] -= 7;
    tree.setNodeValueAtDepth(key, KEY_BIT_WIDTH-1, 1.0);
    EXPECT_EQ(tree.getTreeLogOdds(), 3.0);
    OcTree old_tree (tree);
    tree.setTreeDepth(4);
    EXPECT_EQ(tree.getTreeLogOdds(), 3.0);
    tree.setTreeDepth(KEY_BIT_WIDTH/2);
    EXPECT_EQ(tree.getTreeLogOdds(), 3.0);
    tree.setTreeDepth(KEY_BIT_WIDTH);
    EXPECT_EQ(tree.getTreeLogOdds(), 3.0);
    EXPECT_TRUE(old_tree == tree);
    tree.setTreeDepth(2);
    EXPECT_EQ(tree.getTreeLogOdds(), 3.0);
    tree.setTreeDepth(KEY_BIT_WIDTH);
    tree.clear();
    point3d p2(-1.0, -1.0, -1.0);
    tree.coordToKeyChecked(p2, key);
    tree.setNodeValueAtDepth(key, KEY_BIT_WIDTH-3, 3.0);
    tree.setTreeDepth(1);
    EXPECT_EQ(tree.size(), 0);
  // ------------------------------------------------------------
  } else {
    std::cerr << "Invalid test name specified: " << test_name << std::endl;
    return 1;

  }

  std::cerr << "Test successful.\n";
  return 0;
}
