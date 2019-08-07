/*
 * OctoMap - An Efficient Probabilistic 3D Mapping Framework Based on Octrees
 * http://octomap.github.com/
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


#include <octomap/Abstract3DOcTree.h>
#include <octomap/AbstractOccupancyOcTree.h>
#include <octomap/octomap_types.h>


namespace octomap {

template <class POINT>
Abstract3DOcTree<POINT>::Abstract3DOcTree(double in_resolution)
: OcTreeBaseImpl<OcTreeDataNode<POINT>,AbstractOcTree>(in_resolution)
  {
	//StaticMemberInitializer().ensureLinking();
  }

template <class POINT>
OcTreeDataNode<POINT>* Abstract3DOcTree<POINT>::setNode(const OcTreeKey& key, POINT point, bool lazy_eval) {
	bool createdRoot = false;
	if (this->root == NULL){
		this->root = this->allocNode();
		this->tree_size++;
		createdRoot = true;
	}

	return setNodeRecurs(this->root, createdRoot, key, 0, point, lazy_eval);
}

template <class POINT>
OcTreeDataNode<POINT>* Abstract3DOcTree<POINT>::setNode(const point3d& value, POINT point, bool lazy_eval) {
	OcTreeKey key;
	if (!this->coordToKeyChecked(value, key))
		return NULL;

	return setNode(key, point, lazy_eval);
}

template <class POINT>
OcTreeDataNode<POINT>* Abstract3DOcTree<POINT>::setNode(double x, double y, double z, POINT point, bool lazy_eval) {
	OcTreeKey key;
	if (!this->coordToKeyChecked(x, y, z, key))
		return NULL;

	return setNode(key, point, lazy_eval);
}

// TODO: mostly copy of updateNodeRecurs => merge code or general tree modifier / traversal
template <class POINT>
OcTreeDataNode<POINT>* Abstract3DOcTree<POINT>::setNodeRecurs(OcTreeDataNode<POINT>* node, bool node_just_created, const OcTreeKey& key,
		unsigned int depth, const POINT& point, bool lazy_eval) {
	bool created_node = false;

	assert(node);

	// follow down to last level
	if (depth < this->tree_depth) {
		unsigned int pos = computeChildIdx(key, this->tree_depth -1 - depth);
		if (!this->nodeChildExists(node, pos)) {
			// child does not exist, but maybe it's a pruned node?
			if (!this->nodeHasChildren(node) && !node_just_created ) {
				// current node does not have children AND it is not a new node
				// -> expand pruned node
				this->expandNode(node);
			}
			else {
				// not a pruned node, create requested child
				this->createNodeChild(node, pos);
				created_node = true;
			}
		}

		if (lazy_eval)
			return setNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, point, lazy_eval);
		else {
			OcTreeDataNode<POINT>* retval = setNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, point, lazy_eval);
			// prune node if possible, otherwise set own probability
			// note: combining both did not lead to a speedup!
			if (this->pruneNode(node)){
				// return pointer to current parent (pruned), the just updated node no longer exists
				retval = node;
			}

			return retval;
		}
	}

	// at last level, update node, end of recursion
	else {
		//	      if (use_change_detection) {
		//	        bool occBefore = this->isNodeOccupied(node);
		//	        float valBefore = node->getLogOdds();
		//	        node->setLogOdds(point);
		//
		//	        if (node->getLogOdds() != valBefore || node_just_created)
		//	        {
		//	          valueChangeCallbackWrapper(key, depth, node_just_created, valBefore, occBefore, node->getLogOdds(), this->isNodeOccupied(node));
		//	        }
		//
		//	        if (node_just_created){  // new node
		//	          changed_keys.insert(std::pair<OcTreeKey,bool>(key, true));
		//	        } else if (occBefore != this->isNodeOccupied(node)) {  // occupancy changed, track it
		//	          KeyBoolMap::iterator it = changed_keys.find(key);
		//	          if (it == changed_keys.end())
		//	            changed_keys.insert(std::pair<OcTreeKey,bool>(key, false));
		//	          else if (it->second == false)
		//	            changed_keys.erase(it);
		//	        }
		//	      } else {
		//	        node->setLogOdds(point);
		//	      }
		node->setValue(point);
		return node;
	}
}

}
