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
using namespace std;

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

	num_particles = 512;  // TODO: Set the number of particles
	default_random_engine gen;
	double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta
	std_x = std[0], std_y = std[1], std_theta = std[2]; // Sets standard deviations for x, y, and theta
	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std_x);
	// This line creates a normal (Gaussian) distribution for y
	normal_distribution<double> dist_y(y, std_y);
	// This line creates a normal (Gaussian) distribution for theta
	normal_distribution<double> dist_theta(theta, std_theta);

	for (int i = 0; i < num_particles; ++i)
	{
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
		//<< sample_theta << std::endl;
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
	default_random_engine gen;
	double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta
	std_x = std_pos[0], std_y = std_pos[1], std_theta = std_pos[2]; // Sets standard deviations for x, y, and theta
	// This line creates a normal (Gaussian) distribution for velocity
	//normal_distribution<double> velocity_x_gauss(velocity, std_x);
	// This line creates a normal (Gaussian) distribution for y
	//normal_distribution<double> dist_y(y, std_y);
	// This line creates a normal (Gaussian) distribution for yaw_rate
	//normal_distribution<double> yaw_rate_gauss(yaw_rate, std_theta);
	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(0.0, std_x);
	// This line creates a normal (Gaussian) distribution for y
	normal_distribution<double> dist_y(0.0, std_y);
	// This line creates a normal (Gaussian) distribution for theta
	normal_distribution<double> dist_theta(0.0, std_theta);


	for (int i = 0; i < num_particles; ++i) {
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;
		double current_velocity = velocity; //velocity_x_gauss(gen);
		double current_yaw_rate = yaw_rate; //yaw_rate_gauss(gen);
		if (fabs(current_yaw_rate < 0.000001))
		{
			particles[i].x = p_x + (current_velocity * delta_t * cos(p_theta));
			particles[i].y = p_y + (current_velocity * delta_t * sin(p_theta));
			//particles[i].theta = particles[i].theta;
		}
		else {
			particles[i].x = p_x + ((current_velocity / current_yaw_rate) * (sin(p_theta + (current_yaw_rate * delta_t)) - sin(p_theta)));
			particles[i].y = p_y - ((current_velocity / current_yaw_rate) * (cos(p_theta + (current_yaw_rate * delta_t)) - cos(p_theta)));
			particles[i].theta = particles[i].theta + current_yaw_rate * delta_t;
		}

		// adding random Gaussian noise
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);
	}
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	const std::vector<LandmarkObs> & observations, const Map & map_landmarks) {
	// Updates the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // First, when iterating through each particle, need to transform observation points to map coordinates.
  // Next, associate each observation to its nearest landmark. The distribution can then be calculated.

  // First term of multi-variate normal Gaussian distribution calculated below
  // It stays the same so can be outside the loop
	double std_x_lm, std_y_lm;// Standard deviations for lanmark positions x,y
	std_x_lm = std_landmark[0], std_y_lm = std_landmark[1];// Sets standard deviations for x, y
	double obs_x, obs_y;

	// Iterate through each particle
	for (int i = 0; i < num_particles; ++i) {

		// For calculating multi-variate Gaussian distribution of each observation, for each particle
		double multiv_prob_dist = 1.0;
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		// For each observation
		for (int j = 0; j < observations.size(); ++j) {
			obs_x = observations[j].x;
			obs_y = observations[j].y;
			// transform to map x coordinate
			double x_obs_map;
			x_obs_map = p_x + (cos(p_theta) * obs_x) - (sin(p_theta) * obs_y);

			// transform to map y coordinate
			double y_obs_map;
			y_obs_map = p_y + (sin(p_theta) * obs_x) + (cos(p_theta) * obs_y);

			
			// Find nearest landmark
			vector<Map::single_landmark_s> landmarks = map_landmarks.landmark_list;
			vector<double> landmark_obs_dist(landmarks.size());
			for (int k = 0; k < landmarks.size(); ++k) {

				//looking at those in sensor range of the particle make the valid match
				// If it is in range, put it in the distance vector for calculating nearest neighbor
				double landmark_part_dist = sqrt(pow(p_x - landmarks[k].x_f, 2) + pow(p_y - landmarks[k].y_f, 2));
				if (landmark_part_dist <= sensor_range) {
					landmark_obs_dist[k] = sqrt(pow(x_obs_map - landmarks[k].x_f, 2) + pow(y_obs_map - landmarks[k].y_f, 2));

				}
				else {
					// Need to fill those outside of distance with huge number, or they'll be a zero (and think they are closest)
					landmark_obs_dist[k] = 999999.0;

				}

			}

			// Associate the observation point with its nearest landmark neighbor
			int min_pos = distance(landmark_obs_dist.begin(), min_element(landmark_obs_dist.begin(), landmark_obs_dist.end()));
			float mu_x = landmarks[min_pos].x_f;
			float mu_y = landmarks[min_pos].y_f;

			// Calculate multi-variate Gaussian distribution
			multiv_prob_dist *= multiv_prob(std_x_lm, std_y_lm, x_obs_map, y_obs_map,
				mu_x, mu_y);


		}

		// Update particle weights with combined multi-variate Gaussian distribution
		particles[i].weight = multiv_prob_dist;
		weights[i] = multiv_prob_dist;

	}

}

void ParticleFilter::resample() {
	/**
	 * TODO: Resample particles with replacement with probability proportional
	 *   to their weight.
	 * NOTE: You may find std::discrete_distribution helpful here.
	 *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	 */

	
	//std::discrete_distribution<int> index_val(0, num_particles);

  	//random_device rd;
  	default_random_engine gen;

	vector<Particle> sampled_particles;
	//generate the resampler with spinning weights
	discrete_distribution<int> index(weights.begin(), weights.end());
	for (int i = 0; i < num_particles; i++)
	{
		sampled_particles.push_back(particles[index(gen)]);
	}
	particles = sampled_particles;
}


void ParticleFilter::SetAssociations(Particle & particle,
	const vector<int> & associations,
	const vector<double> & sense_x,
	const vector<double> & sense_y) {
	// particle: the particle to which assign each listed association, 
	//   and association's (x,y) world coordinates mapping
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates
	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;
}


string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
	vector<double> v;

	if (coord == "X") {
		v = best.sense_x;
	}
	else {
		v = best.sense_y;
	}

	std::stringstream ss;
	copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}