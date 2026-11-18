#include "simulation/kinematics.h"
#include <Eigen/Geometry>
#include <iostream>
#include "Eigen/Dense"
#include "acclaim/bone.h"
#include "util/helper.h"
#include <queue>
#include <vector>
#include <cmath>
#include <algorithm>
namespace kinematics {

void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO (FK)
    // You should set these variables:
    //     bone->start_position = Eigen::Vector4d::Zero();
    //     bone->end_position = Eigen::Vector4d::Zero();
    //     bone->rotation = Eigen::Matrix4d::Zero();
    // The sample above just set everything to zero
    // Hint:
    //   1. posture.bone_translations, posture.bone_rotations
    // Note:
    //   1. This function will be called with bone == root bone of the skeleton
    //   2. we use 4D vector to represent 3D vector, so keep the last dimension as "0"
    //   3. util::rotate{Degree | Radian} {XYZ | ZYX}
    //      e.g. rotateDegreeXYZ(x, y, z) means:
    //      x, y and z are presented in degree rotate z degrees along z - axis first, then y degrees along y - axis, and
    //      then x degrees along x - axis
    std::vector<bool> visited(sizeof(posture), 0);
    std::queue<acclaim::Bone*> Q;

    //int init_bone = bone->idx;
    //  dont know root bone's actual end_position(which is also bone1's start_position)
    bone->start_position = posture.bone_translations[bone->idx];
    bone->rotation = util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);  //change to rotation matrix
    bone->end_position = bone->rotation * bone->dir * bone->length + bone->start_position;
    visited[bone->idx] = 1;
    Q.push(bone);
    
    while (!Q.empty()) {
        acclaim::Bone* temp = Q.front();
        visited[temp->idx] = 1;
        //std::cout << temp->idx << std::endl;
        if (temp->idx != 0) {
            // calculate bone idx start_pos v
            temp->start_position = temp->parent->end_position;
            // calculate bone idx rotation matrix
            temp->rotation = temp->parent->rotation * temp->rot_parent_current *
                             util::rotateDegreeZYX(posture.bone_rotations[temp->idx]);
            // calculate bone idx end_pos 
            temp->end_position = temp->rotation * temp->dir * temp->length + temp->start_position;
        }
        
        for (acclaim::Bone* itr = temp->child; itr != nullptr; itr = itr->sibling) {
            if (visited[itr->idx] == 0) {
                visited[itr->idx] == 1;
                Q.push(itr);
                //std::cout << "pushed: " << itr->idx << std::endl;
            }
        }
        Q.pop();
    }
}

Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    // TODO (find x which min(| jacobian * x - target |))
    // Hint:
    //   1. Linear algebra - least squares solution
    //   2. https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Construction
    // Note:
    //   1. SVD or other pseudo-inverse method is useful
    //   2. Some of them have some limitation, if you use that method you should check it.
    Eigen::VectorXd deltatheta(Jacobian.cols());
    deltatheta.setZero();

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Jacobian, Eigen::ComputeThinU | Eigen::ComputeThinV);// (U,Jacobian,V)
    auto jsv = svd.singularValues();
    Eigen::MatrixXd jsvInv(Jacobian.cols(), Jacobian.rows());
    jsvInv.setZero();

    double tolerance = 1e-3; 
    for (int i = 0; i < jsv.size(); ++i) {
        if (jsv(i) < tolerance) jsvInv(i, i) = 0.0;
        else jsvInv(i, i) = 1.0 / jsv(i);
    }

    Eigen::MatrixXd pseudoInverse = svd.matrixV() * jsvInv * svd.matrixU().transpose();
    deltatheta = pseudoInverse * target;
    return deltatheta;
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` (first joint in the chain) will move to.
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param posture The original AMC motion's reference, you need to modify this
 * @param jointChains A 2D vector containing multiple 1D vectors, each of which holds pointers to Eigen::Vector4d
 * constituting a chain.
 * @param boneChains A 2D vector containing multiple 1D vectors, each of which holds pointers to acclaim::Bone
 * constituting a chain.
 * @param currentBasePos The base of the current chain.
 */

bool inverseJacobianIKSolver(std::vector<Eigen::Vector4d> target_pos, acclaim::Bone* end_bone,
                             acclaim::Posture& posture, std::vector<std::vector<Eigen::Vector4d*>>& jointChains,
                             std::vector<std::vector<acclaim::Bone*>>& boneChains, Eigen::Vector4d currentBasePos) {
    constexpr int max_iteration = 1000; //1000
    constexpr double epsilon = 0.1005; //1E-3
    constexpr double step = 0.1; //0.1
    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is the
    // root.
    acclaim::Bone* root_bone = end_bone - end_bone->idx;
    // TODO
    // Perform inverse kinematics (IK)
    // HINTs will tell you what should do in that area.
    // Of course you can ignore it (Any code below this line) and write your own code.
    acclaim::Posture original_posture(posture);

    size_t bone_num = 0;

    // Traverse each chain
    for (int chainIdx = 0; chainIdx < boneChains.size(); ++chainIdx) {
        //std::cout << "it is chainIdx:" << chainIdx << std::endl;
        bone_num = boneChains[chainIdx].size();
        Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);
        Jacobian.setZero();

        for (int iter = 0; iter < max_iteration; ++iter) {
            // std::cout << "it is iter:" << iter << std::endl;
            Eigen::Vector4d desiredVector = target_pos[chainIdx] - *jointChains[chainIdx][0];
            if (desiredVector.norm() < epsilon) {
                break;
            }
            // TODO (compute jacobian)
            //   1. Compute arm vectors
            //   2. Compute jacobian columns, store in `Jacobian`
            // Hint:
            //   1. You should not put rotation in jacobian if it doesn't have that DoF.
            //   2. jacobian.col(/* some column index */) = /* jacobian column */
            int col_idx = 0;
            for (long long i = 0; i < bone_num; i++) {
                Eigen::Vector4d arm_vec = *jointChains[chainIdx][0] - *jointChains[chainIdx][i];
                //std::cout << boneChains[chainIdx][i]->dof << std::endl;
                
                if (boneChains[chainIdx][i]->dofrx) {
                    Jacobian.col(col_idx) = boneChains[chainIdx][i]->rotation.matrix().col(0).cross3(arm_vec);
                    //std::cout << "col_idx:" << col_idx << std::endl;
                    col_idx++;
                }
                if (boneChains[chainIdx][i]->dofry) {
                    Jacobian.col(col_idx) = boneChains[chainIdx][i]->rotation.matrix().col(1).cross3(arm_vec);
                    //std::cout << "col_idx:" << col_idx << std::endl;
                    col_idx++;
                }
                if (boneChains[chainIdx][i]->dofrz) {
                    Jacobian.col(col_idx) = boneChains[chainIdx][i]->rotation.matrix().col(2).cross3(arm_vec);
                    //std::cout << "col_idx:" << col_idx << std::endl;
                    col_idx++;
                }
            }
            //std::cout << "total col_idx:" << col_idx << std::endl;
            //std::cout << "col_idx:" << col_idx << std::endl;

            Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);

            // TODO (update rotation)
            //   Update `posture.bone_rotation` (in euler angle / degrees) using deltaTheta
            // Hint:
            //   1. You can ignore rotation limit of the bone.
            // Bonus:
            //   1. You cannot ignore rotation limit of the bone.
            int delta_index = 0;
            for (long long i = 0; i < bone_num; i++) {
                Eigen::Vector4d current_rotation_deg = posture.bone_rotations[boneChains[chainIdx][i]->idx];
                Eigen::Vector4d delta_rotation;
                delta_rotation.setZero();

                if (boneChains[chainIdx][i]->dofrx) {
                    delta_rotation.x() = deltatheta[delta_index];
                    delta_index++;
                }
                if (boneChains[chainIdx][i]->dofry) {
                    delta_rotation.y() = deltatheta[delta_index];
                    delta_index++;
                }
                if (boneChains[chainIdx][i]->dofrz) {
                    delta_rotation.z() = deltatheta[delta_index];
                    delta_index++;
                }
                Eigen::Vector4d delta_rotation_deg = util::toDegree(delta_rotation);
                Eigen::Vector4d new_rotation_deg = current_rotation_deg + delta_rotation_deg;
                if (boneChains[chainIdx][i]->dofrx) {
                    new_rotation_deg.x() =
                        std::clamp(new_rotation_deg.x(), static_cast<double>(boneChains[chainIdx][i]->rxmin), static_cast<double>(boneChains[chainIdx][i]->rxmax));
                }
                if (boneChains[chainIdx][i]->dofry) {
                    new_rotation_deg.y() =
                        std::clamp(new_rotation_deg.y(), static_cast<double>(boneChains[chainIdx][i]->rymin), static_cast<double>(boneChains[chainIdx][i]->rymax));
                }
                if (boneChains[chainIdx][i]->dofrz) {
                    new_rotation_deg.z() =
                        std::clamp(new_rotation_deg.z(), static_cast<double>(boneChains[chainIdx][i]->rzmin), static_cast<double>(boneChains[chainIdx][i]->rzmax));
                }
                          
                posture.bone_rotations[boneChains[chainIdx][i]->idx] = new_rotation_deg;
                
            }

            forwardSolver(posture, root_bone);
            // Deal with root translation
            if (chainIdx == 0) {
                posture.bone_translations[0] =
                    posture.bone_translations[0] + (currentBasePos - *jointChains[chainIdx][bone_num]);
            }
        }
    }

    // Return whether IK is stable (i.e. whether the ball is reachable) and let the skeleton not swing its hand in the
    // air
    bool stable = true;
    for (int i = 0; i < boneChains.size(); ++i) {
        if ((target_pos[i] - *jointChains[i][0]).norm() > epsilon) {
            stable = false;
        }
    }
    
    // You can replace "!stable" with "false" to see unstable results, but this may lead to some unexpected outcomes.
    if (!stable) {
        posture = original_posture;
        forwardSolver(posture, root_bone);
        return false;
    } else {
        original_posture = posture;
        return true;
    }
}
}  // namespace kinematics
