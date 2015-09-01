


struct particle {
  mass_t mass;
  time_t last_time;
  position_t last_position;
  energy_t last_kinetic_energy;
  energy_t last_energy_of_proximity_to_next;
  velocity_t approx_velocity;
};

void supernaturally_set_velocity(particle& p, velocity_t v) {
  p.approx_velocity = v;
  p.last_kinetic_energy = v*v*p.mass;
}

void adjust_velocity(particle& p) {
  p.approx_velocity = velocity_t(isqrt(p.last_kinetic_energy/p.mass))*sign(p.approx_velocity);
}

const auto energy_of_proximity_constant_spec_units = grams*meters*meters*meters/seconds/seconds;
const auto energy_of_proximity_constant = 10 * energy_of_proximity_constant_spec_units * identity(energy_units*position_units/energy_of_proximity_constant_spec_units);
energy_t energy_of_proximity(position_t const& p0, position_t const& p1) {
  return divide(energy_of_proximity_constant, p1-p0, rounding_strategy<round_down, negative_is_forbidden>());
}
// force_t force_of_proximity(position_t const& p0, position_t const& p1) {
//   return divide(-energy_of_proximity_constant, (p1-p0)*(p1-p0), rounding_strategy<round_down, negative_is_forbidden>());
// }

energy_t invariant_kinetic_energy(particle const& p0, particle const& p1) {
}

struct interaction_info {
  interaction_info(particle const& p0, particle const& p1, time_t new_time) {
    new_pos0 = p0.last_position + p0.approx_velocity*(new_time-p0.last_time);
    
    const energy_t old_eop = p0.last_energy_of_proximity_to_next;
    new_eop = energy_of_proximity(new_pos0, new_pos1);
    
    /*
  with c = s_0m + t_0n
  sm + tn = c
  t = (c-sm)/n
  k = ssm + ttn
  k = ssm + (c-sm)(c-sm)/n
  k = cc/n - 2csm/n + ssmm/n + ssm
  0 = cc/n - 2csm/n + ssmm/n + ssm - k
  0 = ss(m + mm/n) + s(-2cm/n) + (cc/n - k)
  0 = ss(mn + mm) + s(-2cm) + (cc - kn)
  disc = 4(ccmm - (mn + mm)*(cc - kn))
  disc = 4(knmm+knmn-ccmn)
  s = (2cm +- sqrt(disc)) / 2(mn + mm)
  ss = (2cm +- sqrt(disc))(2cm +- sqrt(disc)) / 4(mn + mm)(mn + mm)
  ssm = (2cm +- sqrt(disc))(2cm +- sqrt(disc)) / 4m(n + m)(n + m)
  ssm = (4ccmm +- 4cm*sqrt(disc) + disc) / 4m(n + m)(n + m)
  ssm = (4ccmm +- 4cm*sqrt(disc) + 4(knmm+knmn-ccmn)) / 4m(n + m)(n + m)
  ssm = (ccmm +- cm*sqrt(disc) + (knmm+knmn-ccmn)) / m(n + m)(n + m)
  ssm = (ccm +- c*sqrt(disc) + (knm+knn-ccn)) / (n + m)(n + m)*/
    
    const momentum_t total_momentum = p0.approx_velocity*p0.mass + p1.approx_velocity*p1.mass;
    const  mass_prod = p0.mass*p1.mass;
    //const velocity_t dv = p0.approx_velocity - p1.approx_velocity;
    const  dv_sq = (p0.last_kinetic_energy*p1.mass + p1.last_kinetic_energy*p0.mass)/mass_prod -
      velocity_t(isqrt(p0.last_kinetic_energy*p1.last_kinetic_energy/mass_prod))*sign(p0.approx_velocity)*sign(p1.approx_velocity);
    const  total_momentum_sq = total_momentum*total_momentum;
    const energy_t kinetic_energy_change = old_eop - new_eop;
    const energy_t new_kinetic_energy = p0.last_kinetic_energy + p1.last_kinetic_energy + kinetic_energy_change;
    const  discriminant_quarter = mass_prod*(new_kinetic_energy*(p0.mass+p1.mass) - total_momentum_sq);
    const  discriminant_quarter = mass_prod*(dv_sq*mass_prod + kinetic_energy_change*(p0.mass+p1.mass));
  }
  position_t new_pos0;
  position_t new_pos1;
  energy_t new_eop;
  momentum_t total_momentum;
   total_momentum_sq;
  energy_t new_kinetic_energy;
   discriminant_quarter;
}

void can_interact(particle const& p0, particle const& p1, time_t new_time) {
//   const momentum_t total_momentum = p0.approx_velocity*p0.mass + p1.approx_velocity*p1.mass;
//   const momentum_t total_mass = p0.mass + p1.mass;
//   const momentum_t avg_velocity = total_momentum/total_mass;
//   
  interaction_info info(p0, p1, new_time);
  return (info.discriminant_quarter >= 0);
}

void interact(particle& p0, particle& p1, time_t new_time) {
{
  interaction_info info(p0, p1, new_time);
  assert (info.discriminant_quarter >= 0);
  
  p0.last_time = new_time;
  p0.last_position = info.new_pos0;
  p0.last_kinetic_energy = divide(info.total_momentum_sq*p0.mass - info.total_momentum*isqrt(info.discriminant_quarter*4) + p1.mass*(info.new_kinetic_energy*(p0.mass+p1.mass) - info.total_momentum_sq),
    (p0.mass+p1.mass)*(p0.mass+p1.mass),
    rounding_strategy<round_down, negative_continuous_with_positive>());
  p0.last_energy_of_proximity_to_next = info.new_eop;
  p1.last_time = new_time;
  p1.last_position = info.new_pos1;
  p1.last_kinetic_energy = info.new_kinetic_energy - p0.last_kinetic_energy;
  
  adjust_velocity(p0);
  adjust_velocity(p1);
}

class particles_interact : public event {
public:
  particles_interact(entity_id id0, entity_id id0) : id0(id0),id1(id1) {}
  entity_id id0;
  entity_id id1;

  void operator()(time_steward::accessor* accessor)const override {
    particle& p0 = accessor->get_mut<particle>(accessor->get(id0));
    particle& p1 = accessor->get_mut<particle>(accessor->get(id1));
    interact(p0, p1);
  }
};

class particle_repulsion_preparer : public trigger {
public:
  player_moves_around(entity_id id0, entity_id id0) : id0(id0),id1(id1) {}
  entity_id id0;
  entity_id id1;

  void operator()(time_steward::accessor* accessor)const override {
    particle const& p0 = accessor->get<particle>(accessor->get(id0));
    particle const& p1 = accessor->get<particle>(accessor->get(id1));
    const time_t t = max(p0.last_time, p1.last_time);
    const position_t pos0 = p0.last_position + p0.approx_velocity*(t-p0.last_time);
    const position_t dp = pos1 - pos0;
    const velocity_t dv = p1.approx_velocity - p0.approx_velocity;
    
    time_t t2 = t-1;
    if (dv > 0) {
      t2 = t + divide(dp, dv*10, rounding_strategy<round_down, negative_is_forbidden>());
    }
    if (dv < 0) {
      t2 = t + divide(dp, -dv*12, rounding_strategy<round_down, negative_is_forbidden>());
      while (!can_interact(p0, p1, t2)) {
        t2 = divide(t2+t, 2, rounding_strategy<round_down, negative_continuous_with_positive>());
      }
    }
    
    if (t2 >= t) {
      accessor->anticipate_event(t2, std::shared_ptr<event>(new particles_interact(id0, id1)));
    }
  }
};

class initialize_world : public event {
public:
  void operator()(time_steward::accessor* accessor)const override {
    std::vector<particle> particles;
    const uint32_t num_particles = 4;
    for (uint32_t i = 0; i < num_particles; ++i) {
      particle p;
      p.last_time = 0;
      particles.push_back(p);
    }
    particles[0].mass = 1200 * grams * identity(mass_units/grams);
    particles[0].last_position = 0 * meters * identity(position_units/meters);
    supernaturally_set_velocity(particles[0], 10 * milli*meters/seconds * identity(velocity_units/(milli*meters/seconds)));
    particles[1].mass = 1200 * grams * identity(mass_units/grams);
    particles[1].last_position = 1 * meters * identity(position_units/meters);
    supernaturally_set_velocity(particles[1], -10 * milli*meters/seconds * identity(velocity_units/(milli*meters/seconds)));
    particles[2].mass = 1 * grams * identity(mass_units/grams);
    particles[2].last_position = 2 * meters * identity(position_units/meters);
    supernaturally_set_velocity(particles[2], 0 * milli*meters/seconds * identity(velocity_units/(milli*meters/seconds)));
    particles[3].mass = 8000 * grams * identity(mass_units/grams);
    particles[3].last_position = 20 * meters * identity(position_units/meters);
    supernaturally_set_velocity(particles[3], -10000 * milli*meters/seconds * identity(velocity_units/(milli*meters/seconds)));
    
    for (uint32_t i = 0; i < num_particles; ++i) {
      entity_id id = siphash_id::combining(i);
      entity_ref e = accessor->get(id);
      particle& p = particles[i];
      if (i+1 < num_particles) {
        p.last_energy_of_proximity_to_next = energy_of_proximity(p.last_position, particles[i+1].last_kinetic_energy);
        accessor->set_trigger(id, std::shared_ptr<trigger>(new player_could_hit_walls(id, siphash_id::combining(i+1))));
      }
      accessor->set<particle>(e, p);
    }7
  }
};

