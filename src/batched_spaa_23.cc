/* Code for producing the experiments in Section 6 of 
      "Balanced Allocation in Batches: The Tower of Two Choices"
      by Dimitrios Los and Thomas Sauerwald (SPAA'23)
      [https://arxiv.org/abs/2302.04399]. */
#include <algorithm>
#include <functional>
#include <iostream>
#include <random>

/* Runs the b-Batched setting. This process was introduced in
*     "Multiple-choice balanced allocation in (almost) parallel", 
         by Berenbrink, Czumaj, Englert, Friedetzky, and Nagel (2012)
         [https://arxiv.org/abs/1501.04822].
   
   It starts from an empty load vector and in each round:
     - Allocates b (potentially weighted) balls using the process provided,
       with the load information at the beginning of the batch.
   
   This class keeps track of the load-vector, the maximum load and gap.
   */
template<typename LoadType, typename Generator>
class BatchedSetting {
public:

  /* Initializes the RBB process with the given load vector. */
  BatchedSetting(
        size_t num_bins,
        size_t batch_size,
        const std::function<size_t(const std::vector<LoadType>&, Generator&)> bin_selector,
        const std::function<LoadType(Generator&)> weight_generator)    
    : load_vector_(num_bins, 0), buffer_vector_(num_bins, 0), batch_size_(batch_size), max_load_(0), total_weight_(0),
      bin_selector_(bin_selector), weight_generator_(weight_generator) { 
    
  }
  
  /* Performs an allocation of a batch. */
  void nextRound(Generator& generator) {
    // Phase 1: Perform b allocations.
    size_t n = load_vector_.size();
    
    for (size_t i = 0; i < batch_size_; ++i) {
      size_t idx = bin_selector_(load_vector_, generator);
      LoadType ball_weight = weight_generator_(generator);
      buffer_vector_[idx] += ball_weight;
      total_weight_ += ball_weight;
    }

    // Phase 2: Update and sort the load vector.
    for (int i = 0; i < n; ++i) {
      load_vector_[i] += buffer_vector_[i];
      buffer_vector_[i] = 0;
      max_load_ = std::max(max_load_, load_vector_[i]);
    }
    std::sort(load_vector_.begin(), load_vector_.end(), std::greater<LoadType>());
  }

  /* Returns the current maximum load. */
  LoadType getMaxLoad() const {
    return max_load_;
  }

  /* Returns the current gap. */
  LoadType getGap() const {
    return max_load_ - total_weight_ / LoadType(load_vector_.size());
  }

  /* Returns the current load vector. */
  std::vector<size_t> getLoadVector() const {
    return load_vector_;
  }

private:

  /* Current load vector of the process. */
  std::vector<LoadType> load_vector_;

  /* Buffer vector for the balls allocated in the current batch. */
  std::vector<LoadType> buffer_vector_;

  /* Batch size used in the setting. */
  const size_t batch_size_;

  /* Current maximum load in the load vector. */
  LoadType max_load_;

  /* Total weight of balls in the load vector. */
  LoadType total_weight_;

  /* Function that selects the bin. */
  std::function<size_t(const std::vector<LoadType>&, Generator&)> bin_selector_;

  /* Function that generates the weight of the ball. */
  std::function<LoadType(Generator&)> weight_generator_;
};

/* Sample a bin via the d-Choice process. */
template<class LoadType, class Generator, int d>
size_t d_choice(const std::vector<LoadType>& load_vector, Generator& generator) {
  size_t n = load_vector.size();
  std::uniform_int<size_t> uar(0, n - 1);
  size_t sample = uar(generator);
  for (int i = 1; i < d; ++i) {
    sample = std::max(sample, uar(generator));
  }
  return sample;
}

/* Sample a bin via the Two-Choice process with random tie-breaking. */
template<typename LoadType, typename Generator>
size_t two_choice_with_random_tie_breaking(const std::vector<LoadType>& load_vector, Generator& generator) {
  size_t n = load_vector.size();
  std::uniform_int<size_t> uar(0, n - 1);
  size_t i1 = uar(generator), i2 = uar(generator);
  // Break ties randomly. 
  if (load_vector[i1] <= load_vector[i2]) return i1;
  return i2;
}

/* Sample a bin via the (1+beta)-process. */
template<typename LoadType, typename Generator>
std::function<size_t(const std::vector<LoadType>&, Generator&)> one_plus_beta(double beta) {
  return [beta](const std::vector<LoadType>& load_vector, Generator& generator) {
    size_t n = load_vector.size();
    std::bernoulli_distribution use_two_choices(beta);
    std::uniform_int<size_t> uar(0, n - 1);
    size_t i1 = uar(generator);
    if (!use_two_choices(generator)) {
      return i1;
    }
    size_t i2 = uar(generator);
    return i1 < i2 ? i2 : i1;
  };
}

/* Sample a bin via the (1+beta)-process with random tie-breaking. */
template<typename LoadType, typename Generator>
std::function<size_t(const std::vector<LoadType>&, Generator&)> one_plus_beta_with_random_tie_breaking(double beta) {
  return [beta](const std::vector<LoadType>& load_vector, Generator& generator) {
    size_t n = load_vector.size();
    std::bernoulli_distribution use_two_choices(beta);
    std::uniform_int<size_t> uar(0, n - 1);
    size_t i1 = uar(generator);
    if (!use_two_choices(generator)) {
      return i1;
    }
    size_t i2 = uar(generator);
    if (load_vector[i1] <= load_vector[i2]) return i1;
    return i2;
  };
}

/* Generate a bin selector from the the mixed quantile process,
   which mixes Quantile(delta) process (with probability eta)
   and the One-Choice process. */
template<typename LoadType, typename Generator>
std::function<size_t(const std::vector<LoadType>&, Generator&)> mixed_quantile(double delta, double eta) {
  return [delta, eta](const std::vector<LoadType>& load_vector, Generator& generator) {
    size_t n = load_vector.size();
    std::uniform_int<size_t> uar(0, n - 1);
    std::bernoulli_distribution use_quantile(eta);
    size_t i1 = uar(generator);
    if (!use_quantile(generator) || i1 >= delta * n) {
      return i1;
    }
    return uar(generator);
  };
}

/* Generate unit weight balls. */
template<typename Generator>
size_t unit_weights_generator(Generator& generator) {
  return 1;
}

/* Generate balls with weights from an Exponential distribution. */
template<typename Generator>
double exp_weights_generator(Generator& generator) {
  std::exponential_distribution<double> exp_distribution(1.0);
  return exp_distribution(generator);
}

/* Returns the theoretically optimal parameter for the (1+beta)-process
   in the b-Batched setting. */
double optimal_beta(double beta_factor, int normalised_batch_size, int num_bins) {
  return std::min(1.0, beta_factor * sqrt(1.0 / double(normalised_batch_size) * log(num_bins)));
}

/* Returns the theoretically optimal parameter for the mixed Quantile process. */
double optimal_eta(double eta_factor, int normalised_batch_size, int num_bins) {
  return std::min(1.0, eta_factor * sqrt(1.0 / double(normalised_batch_size) * log(num_bins)));
}

/* Returns the average gap for the given process. */
template<typename LoadType = size_t, typename Generator = std::mt19937>
double avg_gap(
  int num_bins,
  int normalised_batch_size,
  int num_batches,
  int repetitions,
  const std::function<size_t(const std::vector<LoadType>&, Generator&)> process,
  const std::function<LoadType(Generator&)> weights_generator = unit_weights_generator<Generator>) {

  Generator generator;
  double gap_total = 0.0;
  for (int i = 0; i < repetitions; ++i) {
    BatchedSetting<LoadType, std::mt19937> batched_setting(
      num_bins, normalised_batch_size * num_bins, process, weights_generator);
    LoadType cur_gap = 0;
    for (int round = 0; round < num_batches; ++round) {
      batched_setting.nextRound(generator);
      cur_gap = std::max(cur_gap, batched_setting.getGap());
    }
    gap_total += cur_gap;
  }
  return gap_total / double(repetitions);
}

/* Returns the average gap for the One-Choice process. */
double oc_avg_gap(int num_bins, int num_balls, int repetitions) {
  std::mt19937 generator;
  double gap_total = 0.0;
  for (int rep = 0; rep < repetitions; ++rep) {
    std::vector<size_t> load_vector(num_bins, 0);
    std::uniform_int<size_t> uar(0, num_bins - 1);
    for (int j = 0; j < num_balls; ++j) {
      size_t i = uar(generator);
      ++load_vector[i];
    }
    gap_total += *std::max_element(load_vector.begin(), load_vector.end()) - num_balls / double(num_bins);
  }
  return gap_total / double(repetitions);
}

/* Figure 6.1: Average gap for the (1+beta)-process for various beta in (0, 1)
   and batch sizes {20n, 30n, ..., 70n } for m = 20nb.*/
template<typename Generator = std::mt19937>
void generate_gap_vs_beta() {
  int num_bins = 1000;
  int repetitions = 25;
  const std::vector<int> batch_sizes({ 20, 30, 40, 50, 60, 70 });
  for (auto batch_size : batch_sizes) {
    std::cout << "Batch size : " << batch_size << std::endl;
    for (double beta = 0.05; beta <= 1.001; beta += 0.05) {
      std::cout << "(" << beta << ", " << 
        avg_gap(num_bins, batch_size, 20, repetitions, 
          one_plus_beta_with_random_tie_breaking<size_t, Generator>(std::min(beta, 1.0))) << ")" << std::endl;
    }
  }
}

/* Figure 6.2: Average gap for the mixed Quantile process
   with various eta in (0, 1) and batch sizes {20n, 30n, ..., 70n }
   for m = 20nb. */
template<typename Generator = std::mt19937>
void generate_gap_vs_eta() {
  int num_bins = 1'000;
  int repetitions = 25;
  double delta = 0.5;
  const std::vector<int> batch_sizes({ 20, 30, 40, 50, 60, 70 });
  for (auto batch_size : batch_sizes) {
    std::cout << "Batch size : " << batch_size << std::endl;
    for (double eta = 0.05; eta <= 1.001; eta += 0.05) {
      std::cout << "(" << eta << ", " << 
        avg_gap(num_bins, batch_size, num_bins / batch_size, repetitions, 
          mixed_quantile<size_t, Generator>(delta, std::min(eta, 1.0))) << ")" << std::endl;
    }
  }
}

template<typename LoadType, typename Generator=std::mt19937>
void generate_several_processes_vs_batch_size(const std::function<LoadType(Generator&)> weights_generator) {
  int num_bins = 1'000;
  int repetitions = 50;
  std::cout << "Three-Choice:" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 50; ++normalised_batch_size) {
    int num_rounds = num_bins / normalised_batch_size;
    std::cout << "(" << normalised_batch_size << ", " << 
      avg_gap<LoadType, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions,
        d_choice<LoadType, Generator, 3>, weights_generator) << ")" << std::endl;
  }

  std::cout << "Two-Choice:" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 50; ++normalised_batch_size) {
    int num_rounds = num_bins / normalised_batch_size;
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<LoadType, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions,
        two_choice_with_random_tie_breaking<LoadType, Generator>, weights_generator) << ")" << std::endl;
  }

  std::cout << "(1+beta) with beta = 0.5:" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 50; ++normalised_batch_size) {
    int num_rounds = num_bins / normalised_batch_size;
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<LoadType>(num_bins, normalised_batch_size, num_rounds, repetitions,
        one_plus_beta_with_random_tie_breaking<LoadType, Generator>(0.5), weights_generator) << ")" << std::endl;
  }

  std::cout << "(1+beta) with beta = sqrt( (n/b) log n ):" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 50; ++normalised_batch_size) {
    int num_rounds = num_bins / normalised_batch_size;
    double beta = optimal_beta(1.0, normalised_batch_size, num_bins);
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<LoadType>(num_bins, normalised_batch_size, num_rounds, repetitions,
        one_plus_beta_with_random_tie_breaking<LoadType, Generator>(beta), weights_generator) << ")" << std::endl;
  }

  std::cout << "(1+beta) with beta = 0.7 sqrt( (n/b) log n ):" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 50; ++normalised_batch_size) {
    int num_rounds = num_bins / normalised_batch_size;
    double beta = optimal_beta(0.7, normalised_batch_size, num_bins);
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<LoadType>(num_bins, normalised_batch_size, num_rounds, repetitions, 
        one_plus_beta_with_random_tie_breaking<LoadType, Generator>(beta), weights_generator) << ")" << std::endl;
  }
}

/* Figure 6.3: Average gap for various processes for b = {n, 2n, ... , 50n}
   at m = n^2 balls. */
void generate_several_processes_vs_batch_size_unit_weights() {
  generate_several_processes_vs_batch_size<size_t, std::mt19937>(unit_weights_generator<std::mt19937>);
}

/* Figure 6.3: Average gap for various processes for b = {n, 2n, ... , 50n}
   at m = n^2 balls with weights from Exp(1). */
void generate_several_processes_vs_batch_size_exp_weights() {
  generate_several_processes_vs_batch_size<double, std::mt19937>(exp_weights_generator<std::mt19937>);
}

/* Table 6.5: Average gap for batch size b in {20n, 50n, 80n}
   and for n in {10'000, 100'000} at m = 20bn. */
template<typename Generator = std::mt19937>
void generate_table() {
  int repetitions = 20;
  const std::vector<int> ns({ 10'000, 100'000 });
  const std::vector<int> bs({ 20, 50, 80 });
  for (const auto n : ns) {
    std::cout << "Number of bins: " << n << std::endl << std::endl;
    std::cout << "Two-Choice:" << std::endl;
    for (const auto b : bs) {
      std::cout << b << " : " << avg_gap<size_t, Generator>(n, b, 20, repetitions,
        two_choice_with_random_tie_breaking<size_t, Generator>) << std::endl;
    }

    std::cout << "(1+beta) with beta = 0.7 sqrt( (n/b) log n ):" << std::endl;
    for (const auto b : bs) {
      double beta = optimal_beta(0.7, b, n);
      std::cout << b << " : " << avg_gap<size_t, Generator>(n, b, 20, repetitions,
        one_plus_beta_with_random_tie_breaking<size_t, Generator>(beta)) << std::endl;
    }

    std::cout << "Quantile mixed with One-Choice with eta = sqrt( (n/b) log n ):" << std::endl;
    for (const auto b : bs) {
      double eta = optimal_eta(1.0, b, n);
      std::cout << b << " : " << avg_gap<size_t, Generator>(n, b, 20, repetitions,
        mixed_quantile<size_t, Generator>(0.5, eta)) << std::endl;
    }

    std::cout << "One-Choice with m = b balls:" << std::endl;
    for (const auto b : bs) {
      std::cout << b << " : " << oc_avg_gap(n, b * n, repetitions) << std::endl;
    }
  }
}

int main() {
  std::cout << "Figure 6.1:" << std::endl;
  generate_gap_vs_beta();
  std::cout << "Figure 6.2:" << std::endl;
  generate_gap_vs_eta();
  std::cout << "Figure 6.3:" << std::endl;
  generate_several_processes_vs_batch_size_unit_weights();
  std::cout << "Figure 6.4:" << std::endl;
  generate_several_processes_vs_batch_size_exp_weights();
  std::cout << "Table 6.5:" << std::endl;
  generate_table();
  return 0;
}
