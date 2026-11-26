// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Model parameters for the minimal two-orbital Kanamori lattice.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::KanamoriLattice<PointGroup>>> {
public:
  ModelParameters() : t_(1.), U_(4.), U_prime_(2.), J_(0.), J_non_density_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_t() const {
    return t_;
  }
  void set_t(const double t) {
    t_ = t;
  }

  double get_U() const {
    return U_;
  }
  void set_U(const double U) {
    U_ = U;
  }

  double get_U_prime() const {
    return U_prime_;
  }
  void set_U_prime(const double U_prime) {
    U_prime_ = U_prime;
  }

  double get_J() const {
    return J_;
  }
  void set_J(const double J) {
    J_ = J;
  }

  double get_J_non_density() const {
    return J_non_density_;
  }
  void set_J_non_density(const double J) {
    J_non_density_ = J;
  }

private:
  double t_;
  double U_;
  double U_prime_;
  double J_;
  double J_non_density_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::KanamoriLattice<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;
  buffer_size += concurrency.get_buffer_size(t_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(U_prime_);
  buffer_size += concurrency.get_buffer_size(J_);
  buffer_size += concurrency.get_buffer_size(J_non_density_);
  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::KanamoriLattice<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, U_prime_);
  concurrency.pack(buffer, buffer_size, position, J_);
  concurrency.pack(buffer, buffer_size, position, J_non_density_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::KanamoriLattice<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, U_prime_);
  concurrency.unpack(buffer, buffer_size, position, J_);
  concurrency.unpack(buffer, buffer_size, position, J_non_density_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::KanamoriLattice<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  auto tryToRead = [&](const std::string& name, auto& obj) {
    try {
      reader_or_writer.execute(name, obj);
    }
    catch (const std::exception& /*err*/) {
    }
  };

  try {
    reader_or_writer.open_group("Kanamori-model");
    tryToRead("t", t_);
    tryToRead("U", U_);
    tryToRead("U-prime", U_prime_);
    tryToRead("J", J_);
    tryToRead("J-non-density", J_non_density_);
    reader_or_writer.close_group();
  }
  catch (const std::exception& /*r_e*/) {
    std::cout << "\nKanamori-model parameter group is not defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}
