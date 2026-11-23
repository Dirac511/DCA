// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Model parameters for the bilayer two-orbital Kanamori-Hubbard model.
// The hoppings use the same names and meaning as the La3Ni2O7 model; the interaction is the
// Kanamori one built from U and J (U' = U - 2J, same spin U - 3J, exchange J).

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::BilayerKanamori<PointGroup>>> {
public:
  ModelParameters()
      : ex_(0.),
        ez_(0.),
        eb_(0.),
        et_(0.),
        t11x_(-0.4804),
        t11xy_(0.0755),
        t11xx_(-0.0536),
        t22x_(-0.1166),
        t22xy_(-0.0169),
        t22xx_(-0.0145),
        t12x_(0.2367),
        t12xx_(0.),
        s110_(0.00815),
        s11x_(0.),
        s11xy_(0.),
        s11xx_(0.),
        s220_(-0.6296),
        s22x_(0.),
        s22xy_(0.),
        s22xx_(0.),
        s12x_(-0.0257),
        s12xx_(0.),
        U_(4.),
        U_prime_(2.4),
        J_(0.8),
        J_exchange_(0.8) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

#define DCA_BK_GETSET(NAME, MEMBER)  \
  double get_##NAME() const {        \
    return MEMBER;                   \
  }                                  \
  void set_##NAME(const double v) {  \
    MEMBER = v;                      \
  }

  DCA_BK_GETSET(ex, ex_)
  DCA_BK_GETSET(ez, ez_)
  DCA_BK_GETSET(eb, eb_)
  DCA_BK_GETSET(et, et_)
  DCA_BK_GETSET(t11x, t11x_)
  DCA_BK_GETSET(t11xy, t11xy_)
  DCA_BK_GETSET(t11xx, t11xx_)
  DCA_BK_GETSET(t22x, t22x_)
  DCA_BK_GETSET(t22xy, t22xy_)
  DCA_BK_GETSET(t22xx, t22xx_)
  DCA_BK_GETSET(t12x, t12x_)
  DCA_BK_GETSET(t12xx, t12xx_)
  DCA_BK_GETSET(s110, s110_)
  DCA_BK_GETSET(s11x, s11x_)
  DCA_BK_GETSET(s11xy, s11xy_)
  DCA_BK_GETSET(s11xx, s11xx_)
  DCA_BK_GETSET(s220, s220_)
  DCA_BK_GETSET(s22x, s22x_)
  DCA_BK_GETSET(s22xy, s22xy_)
  DCA_BK_GETSET(s22xx, s22xx_)
  DCA_BK_GETSET(s12x, s12x_)
  DCA_BK_GETSET(s12xx, s12xx_)
  DCA_BK_GETSET(U, U_)
  DCA_BK_GETSET(U_prime, U_prime_)
  DCA_BK_GETSET(J, J_)
  DCA_BK_GETSET(J_exchange, J_exchange_)

#undef DCA_BK_GETSET

private:
  double ex_, ez_, eb_, et_;
  double t11x_, t11xy_, t11xx_;
  double t22x_, t22xy_, t22xx_;
  double t12x_, t12xx_;
  double s110_, s11x_, s11xy_, s11xx_;
  double s220_, s22x_, s22xy_, s22xx_;
  double s12x_, s12xx_;
  double U_, U_prime_, J_, J_exchange_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::BilayerKanamori<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;
  for (const double* p :
       {&ex_, &ez_, &eb_, &et_, &t11x_, &t11xy_, &t11xx_, &t22x_, &t22xy_, &t22xx_, &t12x_,
        &t12xx_, &s110_, &s11x_, &s11xy_, &s11xx_, &s220_, &s22x_, &s22xy_, &s22xx_, &s12x_,
        &s12xx_, &U_, &U_prime_, &J_, &J_exchange_})
    buffer_size += concurrency.get_buffer_size(*p);
  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::BilayerKanamori<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  for (const double* p :
       {&ex_, &ez_, &eb_, &et_, &t11x_, &t11xy_, &t11xx_, &t22x_, &t22xy_, &t22xx_, &t12x_,
        &t12xx_, &s110_, &s11x_, &s11xy_, &s11xx_, &s220_, &s22x_, &s22xy_, &s22xx_, &s12x_,
        &s12xx_, &U_, &U_prime_, &J_, &J_exchange_})
    concurrency.pack(buffer, buffer_size, position, *p);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::BilayerKanamori<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  for (double* p :
       {&ex_, &ez_, &eb_, &et_, &t11x_, &t11xy_, &t11xx_, &t22x_, &t22xy_, &t22xx_, &t12x_,
        &t12xx_, &s110_, &s11x_, &s11xy_, &s11xx_, &s220_, &s22x_, &s22xy_, &s22xx_, &s12x_,
        &s12xx_, &U_, &U_prime_, &J_, &J_exchange_})
    concurrency.unpack(buffer, buffer_size, position, *p);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::BilayerKanamori<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  auto tryToRead = [&](const std::string& name, auto& obj) {
    try {
      reader_or_writer.execute(name, obj);
    }
    catch (const std::exception& /*err*/) {
    }
  };

  try {
    reader_or_writer.open_group("BilayerKanamori-model");
    tryToRead("ex", ex_);
    tryToRead("ez", ez_);
    tryToRead("eb", eb_);
    tryToRead("et", et_);
    tryToRead("t11x", t11x_);
    tryToRead("t11xy", t11xy_);
    tryToRead("t11xx", t11xx_);
    tryToRead("t22x", t22x_);
    tryToRead("t22xy", t22xy_);
    tryToRead("t22xx", t22xx_);
    tryToRead("t12x", t12x_);
    tryToRead("t12xx", t12xx_);
    tryToRead("s110", s110_);
    tryToRead("s11x", s11x_);
    tryToRead("s11xy", s11xy_);
    tryToRead("s11xx", s11xx_);
    tryToRead("s220", s220_);
    tryToRead("s22x", s22x_);
    tryToRead("s22xy", s22xy_);
    tryToRead("s22xx", s22xx_);
    tryToRead("s12x", s12x_);
    tryToRead("s12xx", s12xx_);
    tryToRead("U", U_);
    tryToRead("U-prime", U_prime_);
    tryToRead("J", J_);
    tryToRead("J-exchange", J_exchange_);
    reader_or_writer.close_group();
  }
  catch (const std::exception& /*r_e*/) {
    std::cout << "\nBilayerKanamori-model parameter group is not defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}
