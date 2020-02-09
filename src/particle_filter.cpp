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

	num_particles = 128;  // TODO: Set the number of particles
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
		//particles.push_back(Particle());
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1;
		particles.push_back(particle);

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

		if (abs(yaw_rate) == 0)
		{
			particles[i].x = p_x + (velocity * delta_t * cos(p_theta));
			particles[i].y = p_y + (velocity * delta_t * sin(p_theta));
			//particles[i].theta = particles[i].theta;
		}
		else {
			particles[i].x = p_x + ((velocity / yaw_rate) * (sin(p_theta + (yaw_rate * delta_t)) - sin(p_theta)));
			particles[i].y = p_y - ((velocity / yaw_rate) * (cos(p_theta + (yaw_rate * delta_t)) - cos(p_theta)));
			particles[i].theta = particles[i].theta + yaw_rate * delta_t;
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
	std_x_lm = std_landmark[0];
	std_y_lm = std_landmark[1];// Sets standard deviations for x, y


	for (int np = 0; np < num_particles; np++)

	{
		// create a vector for the observations storage.
		vector<LandmarkObs> predicted_observations;
		LandmarkObs obs;

		for (int k = 0; k < observations.size(); k++)

		{

			LandmarkObs trans_obs;

			obs = observations[k];

			trans_obs.x = particles[np].x + obs.x * cos(particles[np].theta) - obs.y * sin(particles[np].theta);

			trans_obs.y = particles[np].y + obs.x * sin(particles[np].theta) + obs.y * cos(particles[np].theta);

			predicted_observations.push_back(trans_obs);

		}
		// make each particles weight to 1.0
		particles[np].weight = 1.0;
		// with this completed the observations conversion with each particle perspective
		// each and every observation from our sensors are handled here.	   		 
		for (int i = 0; i < predicted_observations.size(); i++)

		{

			double valid_dist = sensor_range;

			int assosiated_pos = 0;



			for (int j = 0; j < map_landmarks.landmark_list.size(); j++)

			{

				double nearest_dist = dist(predicted_observations[i].x, predicted_observations[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);

				if (nearest_dist < valid_dist)

				{

					valid_dist = nearest_dist;

					assosiated_pos = j;

				}

			}



			if (assosiated_pos != 0)

			{

				double landmark_obs_x = predicted_observations[i].x;
				double landmark_obs_y = predicted_observations[i].y;
				double nearest_landmark_x = map_landmarks.landmark_list[assosiated_pos].x_f;
				double nearest_landmark_y = map_landmarks.landmark_list[assosiated_pos].y_f;

				//double particle_weight = 1 / (2 * M_PI * std_x_lm * std_y_lm) * exp(-(pow(landmark_obs_x - nearest_landmark_x, 2) / (2 * pow(std_x_lm, 2)) + pow(landmark_obs_y - nearest_landmark_y, 2) / (2 * pow(std_y_lm, 2))));
				double particle_weight = multiv_prob(std_x_lm, std_y_lm, landmark_obs_x, landmark_obs_y, nearest_landmark_x, nearest_landmark_y);
				if (particle_weight > 0)

				{

					particles[np].weight *= particle_weight;

				}

			}

			weights[np] = particles[np].weight;

		}

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