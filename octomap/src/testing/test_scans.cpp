
#include <stdio.h>
#include <octomap/octomap.h>
#include <octomap/math/Utils.h>
#include "testing.h"

using namespace std;
using namespace octomap;

void printUsage(char* self){
  std::cerr << "\nUSAGE: " << self << " spherical_scan.graph  (reference file to compare, required)\n\n";

  exit(1);
}

int main(int argc, char** argv) {
  if (argc != 2){
    printUsage(argv[0]);
  }

  std::string filename = std::string(argv[1]);
  
  ScanGraph referenceGraph;
  EXPECT_TRUE(referenceGraph.readBinary(filename));
  
  // TODO: read in reference graph file


  //##############################################################     

  point3d point_on_surface (4.01, 0.01, 0.01);


  Pointcloud* cloud = new Pointcloud();

  for (int i=-50; i<51; i++) {
    for (int j=-50; j<51; j++) {
      point3d rotated = point_on_surface;
      rotated.rotate_IP(0, DEG2RAD(i*0.5), DEG2RAD(j*0.5));
      cloud->push_back(rotated);
    }
  }
  
  pose6d origin(1.0, 0, -0.5, 0, 0, 0);
  
  ScanGraph graph;
  graph.addNode(cloud, origin); // graph assumes ownership of cloud!
  
  {
    std::cout << "Comparing ScanGraph with reference file at " << filename << std::endl;
    EXPECT_TRUE(graph.size() == referenceGraph.size());
    ScanNode* scanNode = *graph.begin();
    ScanNode* refScanNode = *referenceGraph.begin();
    
    EXPECT_EQ(scanNode->id, refScanNode->id);
    EXPECT_EQ(scanNode->pose, refScanNode->pose);
    EXPECT_EQ(scanNode->scan->size(), refScanNode->scan->size());

    for (size_t i = 0; i < scanNode->scan->size(); ++i){
      EXPECT_FLOAT_EQ((*scanNode->scan)[i].x(), (*refScanNode->scan)[i].x());
      EXPECT_FLOAT_EQ((*scanNode->scan)[i].y(), (*refScanNode->scan)[i].y());
      EXPECT_FLOAT_EQ((*scanNode->scan)[i].z(), (*refScanNode->scan)[i].z());
    }  
    
  }
  // test reading and writing to file
  {
    std::cout << "Testing ScanGraph I/O" << std::endl;
    
    EXPECT_TRUE(graph.writeBinary("spherical_scan_out.graph"));
    
    ScanGraph reReadGraph;
    EXPECT_TRUE(reReadGraph.readBinary("spherical_scan_out.graph"));
    
    EXPECT_TRUE(graph.size() == reReadGraph.size());
    EXPECT_EQ(reReadGraph.size(), 1);
    
    ScanNode* scanNode = *graph.begin();
    ScanNode* readScanNode = *reReadGraph.begin();
    
    EXPECT_EQ(scanNode->id, readScanNode->id);
    EXPECT_EQ(scanNode->pose, readScanNode->pose);
    EXPECT_EQ(scanNode->scan->size(), readScanNode->scan->size());

    for (size_t i = 0; i < scanNode->scan->size(); ++i){
      EXPECT_EQ((*scanNode->scan)[i], (*readScanNode->scan)[i]);
    }
  }
  
  
  // insert into OcTree  
  {
    OcTree tree (0.05);  

    // insert in global coordinates:
    tree.insertPointCloud(*cloud, origin.trans());

    tree.writeBinary("spherical_scan.bt");
  }
  
  cout << "Test done." << endl;
  exit(0);

}
