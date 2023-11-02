#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <vector>
#include <unordered_map>

#include <random>

#include <data_structures.h>


// this will inicate the beginning of the voxel field(x=y=z=0) in world space
extern const float voxel_x_origin;
extern const float voxel_y_origin;
extern const float voxel_z_origin;
// this will adjust voxel size, the voxel size will be voxel_size_scale * 1
extern const float voxel_size_scale;


// set up voxel field
void set_up_voxel_field(voxel_field& V) {
    voxel v1;
    v1.color = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
    v1.density = 50000.0f;
    v1.exist = true;
    for (int i = 2; i < 12; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 6; k++) {
                V.set_voxel(i, j, k, v1);
            }
        }
    }

}

// set up particle system
void set_up_SPH_particles(std::vector<particle>& P) {
    particle p1;
    // p1.currPos = generateRandomVec3();
    p1.prevPos = glm::vec3(0.0f, 0.0f, 0.0f);
    p1.velocity = glm::vec3(0.0f, 0.0f, 0.0f);
    p1.acceleration = glm::vec3(0.0f, 0.0f, 0.0f);
    p1.pamameters = glm::vec3(0.0f, 0.0f, 0.0f);
    p1.deltaCs = glm::vec3(0.0f, 0.0f, 0.0f);

    for (int i = 0; i < P.size(); i++) {
        p1.currPos = generateRandomVec3();
        p1.currPos.y /= 2;
        p1.currPos.y += (y_max - y_min) / 2;
        P[i] = p1;
    }

    //P[0].currPos = glm::vec3(0.125f, 1.0f, 0.125f);
    //P[1].currPos = glm::vec3(0.115f, 2.111f, 0.125f);

}





// definition of the field, contains a 3D array of voxels

voxel_field::voxel_field(int x, int y, int z) {
	x_size = x;
	y_size = y;
	z_size = z;
    NULL_VOXEL.exist = false;
	field.resize(x);
	for (int i = 0; i < x; i++) {
		field[i].resize(y);
		for (int j = 0; j < y; j++) {
			field[i][j].resize(z);
			for (int k = 0; k < z; k++) {
				field[i][j][k].exist = false;
				field[i][j][k].density = 0.0f;
				field[i][j][k].color = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
			}
		}
	}
}
void voxel_field::set_voxel(int x, int y, int z, float density, glm::vec4 color) {
	field[x][y][z].exist = true;
	field[x][y][z].density = density;
	field[x][y][z].color = color;
}
void voxel_field::set_voxel(int x, int y, int z, voxel v) {
	field[x][y][z] = v;
}
voxel & voxel_field::get_voxel(int x, int y, int z) {
    if (x < 0 || x >= x_size || y < 0 || y >= y_size || z < 0 || z >= z_size) {
		return NULL_VOXEL;
	}
	return field[x][y][z];
}
void voxel_field::clear_voxel(int x, int y, int z) {
	field[x][y][z].exist = false;
	field[x][y][z].density = 0.0f;
	field[x][y][z].color = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
}
void voxel_field::clear_all() {
	for (int i = 0; i < x_size; i++) {
		for (int j = 0; j < y_size; j++) {
			for (int k = 0; k < z_size; k++) {
				field[i][j][k].exist = false;
				field[i][j][k].density = 0.0f;
				field[i][j][k].color = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
			}
		}
	}
}
void voxel_field::print_field() {
	for (int i = 0; i < field.size(); i++) {
		for (int j = 0; j < field[i].size(); j++) {
			for (int k = 0; k < field[i][j].size(); k++) {
				if (field[i][j][k].exist) {
					std::cout << "1 ";
				}
				else {
					std::cout << "0 ";
				}
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

    


// avoid some cases that the voxel index is out of bound
void check_voxel_index(int& x, int& y, int& z, voxel_field& V) {
    /*if (x < 0) {
    x = 0;
    }*/
    if (x >= V.x_size) {
        x = V.x_size - 1;
    }
    /*if (y < 0) {
        y = 0;
    }*/
    if (y >= V.y_size) {
        y = V.y_size - 1;
    }
    /*if (z < 0) {
        z = 0;
    }*/
    if (z >= V.z_size) {
        z = V.z_size - 1;
    }
}

// used when need to render the voxel in world space
glm::vec3 voxel_to_world(int x, int y, int z) {
    return glm::vec3(x * voxel_size_scale + voxel_x_origin, y * voxel_size_scale + voxel_y_origin, z * voxel_size_scale + voxel_z_origin);
};
// used in DDA, need to get the x positive face's world pos of the voxel
std::vector<float> voxel_to_world_6_face(int x, int y, int z) {
    std::vector<float> res;// {x_pos,x_neg,y_pos,y_neg,z_pos,z_neg}
    res.push_back(x * voxel_size_scale + voxel_x_origin + voxel_size_scale/2);
    res.push_back(x * voxel_size_scale + voxel_x_origin - voxel_size_scale/2);
    res.push_back(y * voxel_size_scale + voxel_y_origin + voxel_size_scale/2);
    res.push_back(y * voxel_size_scale + voxel_y_origin - voxel_size_scale/2);
    res.push_back(z * voxel_size_scale + voxel_z_origin + voxel_size_scale/2);
    res.push_back(z * voxel_size_scale + voxel_z_origin - voxel_size_scale/2);
    return res;
};
// used when need to get the voxel index from world space(might be a float, which then converted to the int)
std::vector<int> world_to_voxel(glm::vec3 world, voxel_field & V) {
    //std::cout<<"world" << world.x << " " << world.y << " " << world.z << std::endl;
    //glm::vec3 float_val = glm::vec3((world.x - voxel_x_origin) / voxel_size_scale, (world.y - voxel_y_origin) / voxel_size_scale, (world.z - voxel_z_origin) / voxel_size_scale);
    glm::vec3 float_val = glm::vec3((world.x) / voxel_size_scale, (world.y) / voxel_size_scale, (world.z) / voxel_size_scale);
    //std::cout<<"voxel" << float_val.x << " " << float_val.y << " " << float_val.z << std::endl;
    int x = static_cast<int>(std::floor(float_val.x));
    int y = static_cast<int>(std::floor(float_val.y));
    int z = static_cast<int>(std::floor(float_val.z));
    check_voxel_index(x, y, z, V);
    return {x,y,z};
};


// definition of the bounding barrier, contains the max and min coordinates of the box
// particles cannot go beyond the box
bounding_box::bounding_box(GLfloat _x_max, GLfloat _x_min, GLfloat _y_max, GLfloat _y_min, GLfloat _z_max, GLfloat _z_min) {
	this->x_max = _x_max;
	this->x_min = _x_min;
	this->y_max = _y_max;
	this->y_min = _y_min;
	this->z_max = _z_max;
	this->z_min = _z_min;

    // add a small offset because the boundary particles are not exactly 'on' the boundary of the box
    // the true physical boundary is slightly smaller than the rendered box
    GLfloat render_x_max = x_max;
    GLfloat render_x_min = x_min;
    GLfloat render_y_max = y_max;
    GLfloat render_y_min = y_min;
    GLfloat render_z_max = z_max;
    GLfloat render_z_min = z_min;
    render_x_max += 0.11f;
    render_x_min -= 0.11f;
    render_y_max += 0.11f;
    render_y_min -= 0.11f;
    render_z_max += 0.11f;
    render_z_min -= 0.11f;
		
	face_mesh = {
        // Front face
        render_x_min, render_y_min, render_z_max,
        render_x_max, render_y_min, render_z_max,
        render_x_max, render_y_max, render_z_max,
        render_x_min, render_y_min, render_z_max,
        render_x_max, render_y_max, render_z_max,
        render_x_min, render_y_max, render_z_max,
                                        
        // Back face                
        render_x_min, render_y_min, render_z_min,
        render_x_max, render_y_max, render_z_min,
        render_x_max, render_y_min, render_z_min,
        render_x_max, render_y_max, render_z_min,
        render_x_min, render_y_min, render_z_min,
        render_x_min, render_y_max, render_z_min,
                                        
        // Right face               
        render_x_max, render_y_min, render_z_max,
        render_x_max, render_y_min, render_z_min,
        render_x_max, render_y_max, render_z_min,
        render_x_max, render_y_min, render_z_max,
        render_x_max, render_y_max, render_z_min,
        render_x_max, render_y_max, render_z_max,
                                        
        // Left face                
        render_x_min, render_y_min, render_z_min,
        render_x_min, render_y_min, render_z_max,
        render_x_min, render_y_max, render_z_min,
        render_x_min, render_y_max, render_z_min,
        render_x_min, render_y_min, render_z_max,
        render_x_min, render_y_max, render_z_max,
                                        
        // Top face                 
        render_x_min, render_y_max, render_z_max,
        render_x_max, render_y_max, render_z_max,
        render_x_max, render_y_max, render_z_min,
        render_x_min, render_y_max, render_z_max,
        render_x_max, render_y_max, render_z_min,
        render_x_min, render_y_max, render_z_min,

        // Bottom face
        render_x_max, render_y_min, render_z_max,
        render_x_min, render_y_min, render_z_max,
        render_x_max, render_y_min, render_z_min,
        render_x_max, render_y_min, render_z_min,
        render_x_min, render_y_min, render_z_max,
        render_x_min, render_y_min, render_z_min
	};
    GLfloat x_max_edge = render_x_max - 0.005f;
    GLfloat x_min_edge = render_x_min + 0.005f;
    GLfloat y_max_edge = render_y_max - 0.005f;
    GLfloat y_min_edge = render_y_min + 0.005f;
    GLfloat z_max_edge = render_z_max - 0.005f;
    GLfloat z_min_edge = render_z_min + 0.005f;

    face_edge = {
        // Front face
        x_min_edge, y_min_edge, z_max_edge,
        x_max_edge, y_min_edge, z_max_edge,
        x_max_edge, y_min_edge, z_max_edge,
        x_max_edge, y_max_edge, z_max_edge,
        x_max_edge, y_max_edge, z_max_edge,
        x_min_edge, y_max_edge, z_max_edge,
        x_min_edge, y_max_edge, z_max_edge,
        x_min_edge, y_min_edge, z_max_edge,
        // Back face                
        x_min_edge, y_min_edge, z_min_edge,
        x_max_edge, y_min_edge, z_min_edge,
        x_max_edge, y_min_edge, z_min_edge,
        x_max_edge, y_max_edge, z_min_edge,
        x_max_edge, y_max_edge, z_min_edge,
        x_min_edge, y_max_edge, z_min_edge,
        x_min_edge, y_max_edge, z_min_edge,
        x_min_edge, y_min_edge, z_min_edge,
        // Right face                
        x_max_edge, y_min_edge, z_max_edge,
        x_max_edge, y_min_edge, z_min_edge,
        x_max_edge, y_max_edge, z_min_edge,
        x_max_edge, y_max_edge, z_max_edge,
        // Left face                
        x_min_edge, y_min_edge, z_max_edge,
        x_min_edge, y_min_edge, z_min_edge,
        x_min_edge, y_max_edge, z_min_edge,
        x_min_edge, y_max_edge, z_max_edge,

    };

}

// generate a random vec3 in the min and max range
glm::vec3 generateRandomVec3(float _x_max, float _x_min, float _y_max, float _y_min, float _z_max, float _z_min) {
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> distributionX(_x_min, _x_max);
    std::uniform_real_distribution<float> distributionY(_y_min, _y_max);
    std::uniform_real_distribution<float> distributionZ(_z_min, _z_max);

    float randomX = distributionX(rng);
    float randomY = distributionY(rng);
    float randomZ = distributionZ(rng);

    return glm::vec3(randomX * 0.9f, randomY * 0.9f, randomZ * 0.9f);
}






void calculate_SPH_movement(std::vector<particle> & p, float frameTimeDiff, voxel_field & V) {
    int particle_num = p.size();
    // for each particle, calculate the density and pressure
    #pragma omp parallel for
    for (int i = 0; i < particle_num; i++) {
        int cnt = 0;
        float density_sum = 0.f;
        for (int j = 0; j < particle_num; j++) {
            glm::vec3 delta = (p[i].currPos - p[j].currPos);
            float r = length(delta);
            if (r < smoothing_length)
            {
                cnt++;
                density_sum += particle_mass * /* poly6 kernel */ 315.f * glm::pow(smoothing_length * smoothing_length - r * r, 3.f) / (64.f * PI_FLOAT * glm::pow(smoothing_length, 9));
            }
        }
        p[i].pamameters[0] = density_sum;
        p[i].pamameters[1] = glm::max(particle_stiffness * (density_sum - particle_resting_density), 0.f);
        p[i].pamameters[2] = float(cnt);
    }
    // for each particle, calculate the force and acceleration
    #pragma omp parallel for
    for (int i = 0; i < particle_num; i++) {
        glm::vec3 pressure_force = glm::vec3(0.0f, 0.0f, 0.0f);
        glm::vec3 viscosity_force = glm::vec3(0.0f, 0.0f, 0.0f);
        glm::vec3 dCs = glm::vec3(0.0f, 0.0f, 0.0f);
        for (int j = 0; j < particle_num; j++) {
            if (i == j) {
				continue;
			}
            glm::vec3 delta = (p[i].currPos - p[j].currPos);
            float r = length(delta);
            if (r < smoothing_length) {
                if (r == 0.0f) {
					// if the two particles are at the same position, add a small random delta to avoid NaN
                    delta = generateRandomVec3(0.0001f,-0.0001f, 0.0001f, -0.0001f, 0.0001f, -0.0001f);
				}
				// calculate the pressure force
                pressure_force -= particle_mass * (p[i].pamameters[1] + p[j].pamameters[1]) / (2.f * p[j].pamameters[0]) *
                    // gradient of spiky kernel
                    -45.f / (PI_FLOAT * glm::pow(smoothing_length, 6.f)) * glm::pow(smoothing_length - r, 2.f) * glm::normalize(delta);
				// calculate the viscosity force
                viscosity_force += particle_mass * (p[j].velocity - p[i].velocity) / p[j].pamameters[0] *
					// Laplacian of viscosity kernel
					45.f / (PI_FLOAT * glm::pow(smoothing_length, 6.f)) * (smoothing_length - r);
                
                dCs -= particle_mass * glm::pow(smoothing_length * smoothing_length - r * r, 2.f) / p[j].pamameters[0] *
                    // Poly6 kernel
                	945.f / (32.f * PI_FLOAT * glm::pow(smoothing_length, 9.f)) * delta;
			}
            
            


        }


        viscosity_force *= particle_viscosity;
        p[i].acceleration = glm::vec3((pressure_force / p[i].pamameters[0] + viscosity_force / p[i].pamameters[0] + gravity_force));
        p[i].deltaCs = glm::vec3(glm::normalize(dCs));

    }




    // for each particle, calculate the velocity and new position
    #pragma omp parallel for
    for (int i = 0; i < particle_num; i++) {
        glm::vec3 new_velocity = (p[i].velocity + frameTimeDiff * p[i].acceleration);
        glm::vec3 old_velocity = p[i].velocity;
        glm::vec3 old_position = p[i].currPos;
		glm::vec3 new_position = p[i].currPos + frameTimeDiff * new_velocity;
        

   


        p[i].velocity = new_velocity;
        p[i].prevPos = p[i].currPos;
        p[i].currPos = new_position;

        //std::cout<<"old position: "<<old_position.x<<" "<<old_position.y<<" "<<old_position.z<<std::endl;   
        //std::cout<<"new position: "<<new_position.x<<" "<<new_position.y<<" "<<new_position.z<<std::endl;
        
        // ----------------------------------


        // -----------------------particle - voxel collision detection-----------------------
        
        // ------3D-DDA collision------
        //std::cout << "begin DDA" << std::endl;
        // 1. get the initial voxel index & the final voxel index
        std::vector<int> begin_voxel_index = world_to_voxel(p[i].prevPos,V);//check_voxel_index
        std::vector<int> end_voxel_index = world_to_voxel(p[i].currPos,V);
        if (begin_voxel_index == end_voxel_index) {
            voxel& v = V.get_voxel(begin_voxel_index[0], begin_voxel_index[1], begin_voxel_index[2]);
            if (v.exist) {// if inside a voxel, then try to push it out
				// std::cout << "inside collision, bug here" << std::endl;
				// v.color = glm::vec4(1.f, 0.f, 0.f, 1.0f);
                new_velocity.x = -old_velocity.x;
                new_velocity.y = -old_velocity.y;
                new_velocity.z = -old_velocity.z;
                //p[i].currPos = collision_point;
                new_position = old_position;
                new_position += frameTimeDiff * new_velocity;
                p[i].currPos = new_position;

                p[i].velocity = new_velocity;
			}
		}
        else {
            voxel * end_v = &V.get_voxel(end_voxel_index[0], end_voxel_index[1], end_voxel_index[2]);
            voxel * current_v = &V.get_voxel(begin_voxel_index[0], begin_voxel_index[1], begin_voxel_index[2]);
            // DEBUG!
            //current_v->debug = true;
            //end_v->debug = true;

            
            
            std::vector<int> current_voxel_index = begin_voxel_index;
            //std::vector<int> last_voxel_index = current_voxel_index;
            // 2. get the ray direction
            glm::vec3 ray_direction = glm::normalize(p[i].currPos - p[i].prevPos);
            //std::cout << "ray_direction: " << ray_direction.x << " " << ray_direction.y << " " << ray_direction.z << std::endl;
        // 3. get the initial t and final t
            float t_current = 0.f; // begin at 0, start from the beginning of the ray AKA---> the previous position
            float t_end = glm::length(p[i].currPos - p[i].prevPos); // end at the length of the ray
            // 4. get the delta t in each direction X/Y/Z
            float delta_t_x = voxel_size_scale / glm::abs(ray_direction.x);
            float delta_t_y = voxel_size_scale / glm::abs(ray_direction.y);
            float delta_t_z = voxel_size_scale / glm::abs(ray_direction.z);
            int sign_x = ray_direction.x > 0 ? 1 : -1;
            int sign_y = ray_direction.y > 0 ? 1 : -1;
            int sign_z = ray_direction.z > 0 ? 1 : -1;
            // 5. initialize t_next_x, t_next_y, t_next_z
            std::vector<float> voxel_6_face = voxel_to_world_6_face(current_voxel_index[0], current_voxel_index[1], current_voxel_index[2]);
            float t_next_x;
            float t_next_y;
            float t_next_z;
            if (sign_x == 1) {
                t_next_x = glm::abs(voxel_6_face[0] - p[i].prevPos.x) / ray_direction.x;
            }
            else {
                t_next_x = glm::abs(voxel_6_face[1] - p[i].prevPos.x) / ray_direction.x;
            }
            if (sign_y == 1) {
                t_next_y = glm::abs(voxel_6_face[2] - p[i].prevPos.y) / ray_direction.y;
            }
            else {
                t_next_y = glm::abs(voxel_6_face[3] - p[i].prevPos.y) / ray_direction.y;
            }
            if (sign_z == 1) {
                t_next_z = glm::abs(voxel_6_face[4] - p[i].prevPos.z) / ray_direction.z;
            }
            else {
                t_next_z = glm::abs(voxel_6_face[5] - p[i].prevPos.z) / ray_direction.z;
            }


            //std::vector<float> voxel_6_face = voxel_to_world_6_face(current_voxel_index[0], current_voxel_index[1], current_voxel_index[2]);
            //std::cout<<"x_pos"<<voxel_6_face[0]<<voxel_to_world(current_voxel_index[0], current_voxel_index[1], current_voxel_index[2]).x <<std::endl;
            //voxel_to_world_6_face

            if (current_v->exist) {
                //std::cout << "!" << std::endl;
                new_velocity.x = -old_velocity.x;
                new_velocity.y = -old_velocity.y;
                new_velocity.z = -old_velocity.z;
                //p[i].currPos = collision_point;
                new_position = old_position;
                new_position += frameTimeDiff * new_velocity;
                p[i].currPos = new_position;
                p[i].velocity = new_velocity;
            }
            

        // 6. recursively check the voxel in the ray direction, if the voxel is occupied, 
        // then it collides, reverse the velocity and do whatever you want to do
            while (t_current < t_end && !current_v->exist) {
                //std::cout << "DDA" << std::endl;
                //last_voxel_index = current_voxel_index;
                float t_min_next = glm::min(t_next_x, glm::min(t_next_y, t_next_z));
                t_current += t_min_next;
                if (t_min_next == t_next_x)
                {
                    //std::cout << "X" << std::endl;
                    t_next_x += delta_t_x;
                    current_voxel_index[0] += sign_x;
                    if (current_voxel_index[0] < 0 || current_voxel_index[0] >= V.x_size)
                        break;
                }
                else if (t_min_next == t_next_y)
                {
                    //std::cout << "Y" << std::endl;
                    t_next_y += delta_t_y;
                    current_voxel_index[1] += sign_y;
                    if (current_voxel_index[1] < 0 || current_voxel_index[1] >= V.y_size)
                        break;
                }
                else
                {
                    //std::cout << "Z" << std::endl;
                    t_next_z += delta_t_z;
                    current_voxel_index[2] += sign_z;
                    if (current_voxel_index[2] < 0 || current_voxel_index[2] >= V.z_size)
                        break;
                }
                current_v = &V.get_voxel(current_voxel_index[0], current_voxel_index[1], current_voxel_index[2]);
                //current_v->debug = true;
                if (current_v->exist) {
                    //std::cout<<"collide!"<<std::endl;
                    //std::cout << current_voxel_index[0]<<"," << current_voxel_index[1] << "," << current_voxel_index[2] << std::endl;
                    // collide visualization!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    // current_v->color = glm::vec4(1.f, 0.f, 1.f, 1.0f);
                    // get the collision point
                    //glm::vec3 collision_point = p[i].prevPos + (t_current) * ray_direction;
                    // std::cout << "collision point" << collision_point.x << " " << collision_point.y << " " << collision_point.z << std::endl;

                    new_velocity.x = -old_velocity.x;
                    new_velocity.y = -old_velocity.y;
                    new_velocity.z = -old_velocity.z;
                    //new_position = collision_point;
                    new_position = old_position;
                    p[i].currPos = new_position;

                    p[i].velocity = new_velocity*0.8f;
                    break;
                }
            }
            //std::cout << "after c-d position: " << p[i].currPos.x << " " << p[i].currPos.y << " " << p[i].currPos.z << std::endl;
            





            

        }

        // check collision with the bounding box
        if (new_position.y < y_min)
        {
            new_position.y = y_min;
            new_velocity.y *= -1 * wall_damping;
        }
        else if (new_position.y > y_max)
        {
            new_position.y = y_max;
            new_velocity.y *= -1 * wall_damping;
        }
        if (new_position.x < x_min)
        {
            new_position.x = x_min;
            new_velocity.x *= -1 * wall_damping;
        }
        else if (new_position.x > x_max)
        {
            new_position.x = x_max;
            new_velocity.x *= -1 * wall_damping;
        }
        if (new_position.z < z_min)
        {
            new_position.z = z_min;
            new_velocity.z *= -1 * wall_damping;
        }
        else if (new_position.z > z_max)
        {
            new_position.z = z_max;
            new_velocity.z *= -1 * wall_damping;
        }

        p[i].velocity = new_velocity;
        p[i].currPos = new_position;

        //std::cout << "final position: " << p[i].currPos.x << " " << p[i].currPos.y << " " << p[i].currPos.z << std::endl;




        // // ------simplest collision detection, just reverse the velocity if this pos has a voxel------
        // std::vector<int> voxel_index = world_to_voxel(new_position,V);
        // int x = voxel_index[0];
        // int y = voxel_index[1];
        // int z = voxel_index[2];
        // // avoid some cases that the voxel index is out of bound
        // /*if (x < 0) {
		// 	x = 0;
		// }*/
        // 
        // //std::cout << x << " " << y << " " << z << std::endl;
        // voxel & v = V.get_voxel(x, y, z);
        // if (v.exist) {
        //     
        //     std::cout << "collision" << std::endl;
        //     v.color = glm::vec4(1.f, 0.f, 0.f, 1.0f);
        //     new_velocity.x = -new_velocity.x;
        //     new_velocity.y = -new_velocity.y;
        //     new_velocity.z = -new_velocity.z;
        //     p[i].velocity = new_velocity;
        // 
        // }

    }
}

void calculate_voxel_erosion(std::vector<particle>& p, float frameTimeDiff, voxel_field& V) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < V.x_size; i++) {
        for (int j = 0; j < V.y_size; j++) {
            for (int k = 0; k < V.z_size; k++) {
                voxel * v = &V.get_voxel(i,j,k);
                if (v->exist) {
                    for (int n = 0; n < particle_num; n++) {
                        glm::vec3 delta = (p[n].currPos - voxel_to_world(i, j, k));
						float r = length(delta);
                        if (r < smoothing_length) {
							v->density -= particle_mass * /* poly6 kernel */ 315.f * glm::pow(smoothing_length * smoothing_length - r * r, 3.f) / (64.f * PI_FLOAT * glm::pow(smoothing_length, 9));
						} 
                        if (v->density < voxel_density_threshold) {
                            v->exist = false;
                        }

                    }


					//v->color = glm::vec4(0.f, 0.f, 0.f, 1.0f);
				}

            }
        }
    }



}