/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
	if (is_initialized) {

		return;

	}

	num_particles = 1000;  // TODO: Set the number of particles
	std::default_random_engine gen;
	double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta
	std_x = std[0], std_y = std[1], std_theta = std[2]; // Sets standard deviations for x, y, and theta
	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std_x);
	// This line creates a normal (Gaussian) distribution for y
	normal_distribution<double> dist_y(y, std_y);
	// This line creates a normal (Gaussian) distribution for theta
	normal_distribution<double> dist_theta(theta, std_theta);

	
	for (int i = 0; i < num_particles; ++i) {
		double sample_x, sample_y, sample_theta;

		// TODO: Sample from these normal distributions like this: 
		//   sample_x = dist_x(gen);
		//   where "gen" is the random engine initialized earlier.
		particles.push_back(Particle());
		particles[i].id = i;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
		particles[i].weight = 1;

		// Print your samples to the terminal.
		//std::cout << "Sample " << i + 1 << " " << sample_x << " " << sample_y << " "
		//	<< sample_theta << std::endl;
		weights.push_back(1);//set weights of all particles to 1

	}
	
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
	std::default_random_engine gen;
	double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta
	std_x = std_pos[0], std_y = std_pos[1], std_theta = std_pos[2]; // Sets standard deviations for x, y, and theta
	// This line creates a normal (Gaussian) distribution for velocity
	normal_distribution<double> velocity_x_gauss(velocity, std_x);
	// This line creates a normal (Gaussian) distribution for y
	//normal_distribution<double> dist_y(y, std_y);
	// This line creates a normal (Gaussian) distribution for yaw_rate
	normal_distribution<double> yaw_rate_gauss(yaw_rate, std_theta);


	for (int i = 0; i < num_particles; ++i) {
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;
		double current_velocity = velocity_x_gauss(gen);
		double current_yaw_rate = yaw_rate_gauss(gen);
		if (current_yaw_rate == 0.0)
		{
			particles[i].x = p_x + (current_velocity * delta_t * cos(p_theta));
			particles[i].y = p_y + (current_velocity * delta_t * sin(p_theta));
			particles[i].theta = particles[i].theta;
		}
		else {
			particles[i].x = p_x + current_velocity * (sin(p_theta+(current_yaw_rate *delta_t))-sin(current_yaw_rate)) / current_yaw_rate;
			particles[i].y = p_y - current_velocity * (cos(p_theta + (current_yaw_rate * delta_t)) - cos(current_yaw_rate)) / current_yaw_rate;
			particles[i].theta = particles[i].theta + current_yaw_rate * delta_t;
		}

	}



}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
   // sort the vectors
	//use dist function
	// known pitfall is might it may locate same observation for two predcitions.
	double nearest_dist;
	double save_id=0;
	double past_dist = 65535.0;

	for (int i = 0; i <= predicted.size(); i++) {
		for (int j = 0; j <= observations.size(); j++)
		{
			nearest_dist = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
			if (nearest_dist < past_dist)
			{
				save_id = j;
				past_dist = nearest_dist;
			}
			
		}
		predicted[i].id = save_id;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
	std::default_random_engine gen;
	double std_x_lm, std_y_lm;// Standard deviations for lanmark positions x,y
	std_x_lm = std_landmark[0], std_y_lm = std_landmark[1];// Sets standard deviations for x, y
	double obs_x, obs_y;
	double noisy_obs_x, noisy_obs_y;
	vector<LandmarkObs> predicted_observations;
	vector<LandmarkObs> map_observations;
	double particle_weight = 1;
	for (int i = 0; i < num_particles; ++i) {
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		for (int j = 0; j <= observations.size(); j++) {

			obs_x = observations[j].x;
			obs_y = observations[j].y;
			// This line creates a normal (Gaussian) distribution for x
			//normal_distribution<double> dist_x(obs_x, std_x_lm);
			// This line creates a normal (Gaussian) distribution for y
			//normal_distribution<double> dist_y(obs_y, std_y_lm);

			noisy_obs_x = obs_x; // dist_x(gen);
			noisy_obs_y = obs_y; // dist_y(gen);

			// transform to map x coordinate
			double x_obs_map;
			x_obs_map = p_x + (cos(p_theta) * noisy_obs_x) - (sin(p_theta) * noisy_obs_y);

			// transform to map y coordinate
			double y_obs_map;
			y_obs_map = p_y + (sin(p_theta) * noisy_obs_x) + (cos(p_theta) * noisy_obs_y);

			LandmarkObs conv_obs;
			conv_obs.x = x_obs_map;
			conv_obs.y = y_obs_map;
			predicted_observations.push_back(conv_obs);
		}
		for (int k = 0; k < map_landmarks.landmark_list.size(); k++)
		{
			double map_lm_x = map_landmarks.landmark_list[k].x_f;
			double map_lm_y = map_landmarks.landmark_list[k].y_f;
			int map_lm_id = map_landmarks.landmark_list[k].id_i;
			LandmarkObs map_lm;
			if ((map_lm_x - p_x <= sensor_range) && (map_lm_y - p_y <= sensor_range))
			{
				map_lm.id = map_lm_id;
				map_lm.x = map_lm_x;
				map_lm.y = map_lm_y;
				map_observations.push_back(map_lm);
			}


		}
		ParticleFilter::dataAssociation(predicted_observations, map_observations);

		for (int mn = 0; mn <= predicted_observations.size(); mn++)
		{
			int id = predicted_observations[mn].id;
			if (id > 0)
			{
				particle_weight *= multiv_prob(std_x_lm, std_y_lm, predicted_observations[mn].x, predicted_observations[mn].y,
					map_observations[id].x, map_observations[id].y);
			}
			
		}

		particles[i].weight = particle_weight;
		weights[i] = particle_weight;
	}
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

	std::default_random_engine gen;
	//std::discrete_distribution<int> index_val(0, num_particles);
	
	double beta = 0.0;
	int index;
	double max_w = *max_element(weights.begin(), weights.end());

	std::uniform_real_distribution<double> beta_val(0.0, 2*max_w);
	std::uniform_int_distribution<int> distIndex(0, num_particles - 1);
	
	vector<Particle> sampled_particles;
	//generate the resampler with spiining weights
	for (int i = 0; i < num_particles; i++)
	{
		index = distIndex(gen);
		beta += beta_val(gen);
		while (beta > weights[index])
		{
			beta -= weights[index];
			index = (index + 1) % num_particles;

		}
		sampled_particles.push_back(particles[index]);
	}
	particles = sampled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}