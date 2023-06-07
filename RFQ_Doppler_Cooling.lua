simion.workbench_program()

-- Variables adjustable during flight:

adjustable _amu_mass_per_charge   =   100.0  -- mass/charge tune point (u/e)
                                             -- (particles of this m/z pass)
adjustable pe_update_each_usec    = 0.05   -- potential energy display
                                             -- update period (microsec)
                                             -- (for display purposes only)

-- Variables adjustable only at beginning of flight:

adjustable effective_radius_in_cm   = 0.40   -- half the minimum distance between
                                             -- opposite rods (cm)
adjustable phase_angle_deg          = 0.0    -- quad entry phase angle of ion (deg)
adjustable frequency_hz             = 1.1E6  -- RF frequency of quad (Hz)

adjustable q_u = 0.620
adjustable a_u = 0.122

--------------------------------------------------------------------------------

file = io.open("SR_v_long.txt", "w")
file_x = io.open("SR_v_long_x.txt", "w")
file_y = io.open("SR_v_long_y.txt", "w")
file_z = io.open("SR_v_long_z.txt", "w")

adjustable T1_min = 23
adjustable T1_max = 23
adjustable T2_min = 34
adjustable T2_max = 34

adjustable V_endcap_min = 0
adjustable V_endcap_max = 0 --30
adjustable V_add = 4 --4


local x_0 = 0
local y_0 = 0
local z_0 = 0
local x = 20
local y = 20
local z = 46.5

local c = 299792.458 -- mm/mus
local omega_resonance = 710408668246445.5 -- Hz
local Gamma = 135*10^6 -- Hz 

local omega_laser_x_list = {710403861442020.1}
local omega_laser_y_list = {710403861442020.1}
local omega_laser_z_list = {710405046262198.2}

local omega_laser_x = omega_laser_x_list[1]
local omega_laser_y = omega_laser_y_list[1]
local omega_laser_z = omega_laser_z_list[1]

local I_laser = 2*37.2 -- mW/cm^2
local I_sat = 37.2 -- mW/cm^2
local lamda = 4.22*10^(-4) -- mm
local lifetime = 7.39*10^(-3) -- mus
local Gamma = 135*10^6 -- Hz 
--local photon_mom = 6.626*6.022*10^(-4)/4.22 -- mm*amu/mus
local photon_mom = 9.4569962697883*10^-4 -- mm*amu/mus

local i = 0
local j = 0
local k = 0
local lol = 0

adjustable i_clone = 0
adjustable j_clone = 0
adjustable k_clone = 0

local counter = 0

function segment.flym() -- Called at the beginning of every flym
  sim_trajectory_image_control = 1 -- Don't preserve trajectories


  for i = T1_min, T1_max, 0.1 do
    for j = T2_min, T2_max, 0.1 do
      for lol = 0, 0, 1 do
        i_clone = i
        j_clone = j
        k_clone = k
        run()
      end 
    end
  end
end


-- Temporary variables used internally.
local scaled_rf  -- a factor used in the RF component
local omega      -- frequency_hz (reexpressed in units of radians/usec)
local theta      -- phase_angle_deg (reexpressed in units of radians)
local last_pe_update = 0.0 -- last potential energy surface update time (usec)

-- SIMION segment called by SIMION to set adjustable electrode voltages
-- in the current potential array instance.
-- NOTE: this is called frequently, multiple times per time-step (by
-- Runge-Kutta), so performance concerns here can be important.



function segment.fast_adjust()
  -- See "Overview of Quad Equations" comments for details.

  --print(count,ion_time_step, ion_time_of_flight)
  if not scaled_rf then
      -- Initialize constants if not already initialized.
      -- These constants don't change during particle flight,
      -- so we can calculate them once and reuse them.
      -- Reusing them is a bit more efficient (~25% by one estimate)
      -- than recalculating them on every fast_adjust call.
      scaled_rf = effective_radius_in_cm^2 * frequency_hz^2 * 1.022442E-11 * q_u
      theta = phase_angle_deg * (math.pi / 180)
      omega = frequency_hz * (1E-6 * 2 * math.pi)
  end

  --print(i,j)

  local rfvolts = 2*scaled_rf * _amu_mass_per_charge -- *2
  --print("The V_RF is = " .. rfvolts)
  local dcvolts = 0.5*rfvolts * (a_u/(2*q_u) * 0.5*2)
  --print("The V_DC is = " .. dcvolts)
  local tempvolts = sin(ion_time_of_flight * omega + theta) * rfvolts + dcvolts
  --print("The V_total = " .. tempvolts)

  -- Finally, apply adjustable voltages to rod electrodes.

  if (ion_time_of_flight <= j_clone) then

    adj_elect01 = tempvolts + 0.25*5.7  
    adj_elect02 = - dcvolts + 0.25*5.7 

    adj_elect03 = tempvolts + 0.5*5.7 
    adj_elect04 = - dcvolts + 0.5*5.7  

    adj_elect05 = tempvolts + 0.75*5.7 
    adj_elect06 = - dcvolts + 0.75*5.7 

    adj_elect07 = tempvolts + 1.0*5.7  
    adj_elect08 = - dcvolts + 1.0*5.7  

    adj_elect09 = tempvolts + 1.25*5.7  
    adj_elect10 = - dcvolts + 1.25*5.7  

    adj_elect11 = tempvolts + 1.50*5.7  
    adj_elect12 = - dcvolts + 1.50*5.7  

    adj_elect13 = tempvolts + 1.75*5.7  
    adj_elect14 = - dcvolts + 1.75*5.7  

    adj_elect15 = tempvolts + 2.0*5.7  
    adj_elect16 = - dcvolts + 2.0*5.7

    adj_elect17 = tempvolts + 2.0*5.7  
    adj_elect18 = - dcvolts + 2.0*5.7

    adj_elect19 = tempvolts + 2.0*5.7 
    adj_elect20 = - dcvolts + 2.0*5.7

  elseif (ion_time_of_flight > j_clone) then

    adj_elect01 = tempvolts 
    adj_elect02 = - dcvolts

    adj_elect03 = tempvolts
    adj_elect04 = - dcvolts 

    adj_elect05 = tempvolts
    adj_elect06 = - dcvolts

    adj_elect07 = tempvolts
    adj_elect08 = - dcvolts

    adj_elect09 = tempvolts  
    adj_elect10 = - dcvolts 

    adj_elect11 = tempvolts  
    adj_elect12 = - dcvolts  

    adj_elect13 = tempvolts 
    adj_elect14 = - dcvolts

    adj_elect15 = tempvolts + V_add
    adj_elect16 = - dcvolts + V_add

    adj_elect17 = tempvolts  
    adj_elect18 = - dcvolts

    adj_elect19 = tempvolts + V_add
    adj_elect20 = - dcvolts + V_add
  end

end


--local omega_laser_z = laser_list[1]

--------------------------------------------------------------------------------

-- SIMION segment called by SIMION after every time-step.
function segment.other_actions()

  local KE = 1.05 * 0.5 * (ion_vx_mm^2 + ion_vy_mm^2 + ion_vz_mm^2)
  --local KE = speed_to_ke(ion_vz_mm,100)

  -- Update potential energy surface display periodically.
  -- The performance overhead of this in non-PE views is only a few percent.
  -- NOTE: the value inside abs(...) can be negative when a new ion is flown.
  if abs(ion_time_of_flight - last_pe_update) >= pe_update_each_usec then
    last_pe_update = ion_time_of_flight
    sim_update_pe_surface = 1    -- Request a PE surface display update.
  end

  if (ion_splat ~= 0 ) then -- Unstable trajectory
    --print("UNSTABLE")
    --mark()
    counter = 0
  elseif (ion_time_of_flight > 1000*10^3) then
    --print("STABLE")
    --mark()
    ion_splat = - 4 -- Ion killed
    counter = 0
  end

  rand1 = rand()
  rand2 = rand()
  gaus_rand = math.sqrt(-2*math.log(rand1))*math.cos(2*math.pi*rand2) -- Mean: 0 STDEV: 1

  if (ion_time_of_flight < 295*10^3*0.1) then
    omega_laser_x = omega_laser_x_list[1] --+ gaus_rand * 1 * 10^6
    omega_laser_y = omega_laser_y_list[1] --+ gaus_rand * 1 * 10^6
    omega_laser_z = omega_laser_z_list[1] --+ gaus_rand * 1 * 10^6
  -- elseif (ion_time_of_flight > 295*10^3*0.1 and ion_time_of_flight < 295*10^3*0.2) then
  --   omega_laser_x = omega_laser_x_list[2] + gaus_rand * 1 * 10^6
  --   omega_laser_y = omega_laser_y_list[2] + gaus_rand * 1 * 10^6
  --   omega_laser_z = omega_laser_z_list[2] + gaus_rand * 1 * 10^6
  -- elseif (ion_time_of_flight > 295*10^3*0.2 and ion_time_of_flight < 295*10^3*0.3) then
  --   omega_laser_x = omega_laser_x_list[3] + gaus_rand * 1 * 10^6
  --   omega_laser_y = omega_laser_y_list[3] + gaus_rand * 1 * 10^6
  --   omega_laser_z = omega_laser_z_list[3] + gaus_rand * 1 * 10^6
  -- elseif (ion_time_of_flight > 295*10^3*0.3 and ion_time_of_flight < 295*10^3*0.4) then
  --   omega_laser_x = omega_laser_x_list[4] + gaus_rand * 1 * 10^6
  --   omega_laser_y = omega_laser_y_list[4] + gaus_rand * 1 * 10^6
  --   omega_laser_z = omega_laser_z_list[4] + gaus_rand * 1 * 10^6
  -- elseif (ion_time_of_flight > 295*10^3*0.4 and ion_time_of_flight < 295*10^3*0.5) then
  --   omega_laser_x = omega_laser_x_list[5] + gaus_rand * 1 * 10^6
  --   omega_laser_y = omega_laser_y_list[5] + gaus_rand * 1 * 10^6
  --   omega_laser_z = omega_laser_z_list[5] + gaus_rand * 1 * 10^6
  -- elseif (ion_time_of_flight > 295*10^3*0.5 and ion_time_of_flight < 295*10^3*0.6) then
  --   omega_laser_x = omega_laser_x_list[6] + gaus_rand * 1 * 10^6
  --   omega_laser_y = omega_laser_y_list[6] + gaus_rand * 1 * 10^6
  --   omega_laser_z = omega_laser_z_list[6] + gaus_rand * 1 * 10^6
  -- elseif (ion_time_of_flight > 295*10^3*0.6 and ion_time_of_flight < 295*10^3*0.7) then
  --   omega_laser_x = omega_laser_x_list[7] + gaus_rand * 1 * 10^6
  --   omega_laser_y = omega_laser_y_list[7] + gaus_rand * 1 * 10^6
  --   omega_laser_z = omega_laser_z_list[7] + gaus_rand * 1 * 10^6
  -- elseif (ion_time_of_flight > 295*10^3*0.7 and ion_time_of_flight < 295*10^3*0.8) then
  --   omega_laser_x = omega_laser_x_list[8] + gaus_rand * 1 * 10^6
  --   omega_laser_y = omega_laser_y_list[8] + gaus_rand * 1 * 10^6
  --   omega_laser_z = omega_laser_z_list[8] + gaus_rand * 1 * 10^6
  -- elseif (ion_time_of_flight > 295*10^3*0.8 and ion_time_of_flight < 295*10^3*0.9) then
  --   omega_laser_x = omega_laser_x_list[9] + gaus_rand * 1 * 10^6
  --   omega_laser_y = omega_laser_y_list[9] + gaus_rand * 1 * 10^6
  --   omega_laser_z = omega_laser_z_list[9] + gaus_rand * 1 * 10^6
  -- elseif (ion_time_of_flight > 295*10^3*0.9 and ion_time_of_flight < 295*10^3*1.0) then
  --   omega_laser_x = omega_laser_x_list[10] + gaus_rand * 1 * 10^6
  --   omega_laser_y = omega_laser_y_list[10] + gaus_rand * 1 * 10^6
  --   omega_laser_z = omega_laser_z_list[10] + gaus_rand * 1 * 10^6
  end


  -- Doppler Laser Cooling

  if (ion_time_of_flight > j_clone + 20 + 100) then

    local delta_x = omega_laser_x - omega_resonance - ion_vx_mm * omega_laser_x / c 
    local delta_y = omega_laser_y - omega_resonance - ion_vy_mm * omega_laser_y / c 
    local delta_z = omega_laser_z - omega_resonance - ion_vz_mm * omega_laser_z / c 

    local R_x = 0.5*Gamma*(I_laser/I_sat)/(I_laser/I_sat + 1 + (2*delta_x/Gamma)^2) -- Hz
    local R_x_per_dt = R_x*7.39/(10^9) / 3
    local R_y = 0.5*Gamma*(I_laser/I_sat)/(I_laser/I_sat + 1 + (2*delta_y/Gamma)^2) -- Hz
    local R_y_per_dt = R_y*7.39/(10^9) / 3
    local R_z = 0.5*Gamma*(I_laser/I_sat)/(I_laser/I_sat + 1 + (2*delta_z/Gamma)^2) -- Hz
    local R_z_per_dt = R_z*7.39/(10^9) / 3

    function sign(p_0,p,ion_v) -- Should also be valid for the 3D case.
      if ((p-p_0)>0 and ion_v<0) then
        return 1
      end
      if ((p-p_0)>0 and ion_v>0) then
        return 1
      end
      if ((p-p_0)<0 and ion_v<0) then
        return -1
      end
      if ((p-p_0)<0 and ion_v>0) then
        return -1
      end
      if ((p-p_0)>0 and ion_v==0) then -- Extra cases
        return 1
      end
      if ((p-p_0)<0 and ion_v==0) then
        return -1
      end
      if ((p-p_0)==0) then 
        return 0 
      end
    end 

    counter = counter + 1

    if counter == 1000 then
      -- file:write(ion_time_of_flight, ',') 
      -- file:write(KE, ',') 
      -- file:write(R_x_per_dt, ',') 
      -- file:write(R_y_per_dt, ',')
      -- file:write(R_z_per_dt, '\n')
      mark()
      counter = 0
    end
      
    local random_number = rand()

    if (R_x_per_dt > random_number) then

      ion_vx_mm = ion_vx_mm + sign(x_0,x,ion_vx_mm) * photon_mom / _amu_mass_per_charge

      local KE_x = speed_to_ke(ion_vx_mm,100)
      -- file_x:write(ion_time_of_flight, ',') 
      -- file_x:write(KE_x, ',') 
      -- file_x:write(R_x_per_dt, '\n')

      local random_number_emit = rand()

      if (random_number_emit <= 1/6 and random_number_emit > 0) then
        ion_vx_mm = ion_vx_mm + photon_mom / _amu_mass_per_charge
      elseif (random_number_emit <= 2/6 and random_number_emit > 1/6) then
        ion_vx_mm = ion_vx_mm - photon_mom / _amu_mass_per_charge
      elseif (random_number_emit <= 3/6 and random_number_emit > 2/6) then
        ion_vy_mm = ion_vy_mm + photon_mom / _amu_mass_per_charge
      elseif (random_number_emit <= 4/6 and random_number_emit > 3/6) then
        ion_vy_mm = ion_vy_mm - photon_mom / _amu_mass_per_charge
      elseif (random_number_emit <= 5/6 and random_number_emit > 4/6) then
        ion_vz_mm = ion_vz_mm + photon_mom / _amu_mass_per_charge
      elseif (random_number_emit <= 1 and random_number_emit > 5/6) then
        ion_vz_mm = ion_vz_mm - photon_mom / _amu_mass_per_charge
      end


    end

    local random_number_1 = rand()

    if (R_y_per_dt > random_number_1) then

      ion_vy_mm = ion_vy_mm + sign(y_0,y,ion_vy_mm) * photon_mom / _amu_mass_per_charge

      local KE_y = speed_to_ke(ion_vy_mm,100)
      -- file_y:write(ion_time_of_flight, ',') 
      -- file_y:write(KE_y, ',') 
      -- file_y:write(R_y_per_dt, '\n')

      local random_number_emit_1 = rand()

      if (random_number_emit_1 <= 1/6 and random_number_emit_1 > 0) then
        ion_vx_mm = ion_vx_mm + photon_mom / _amu_mass_per_charge
      elseif (random_number_emit_1 <= 2/6 and random_number_emit_1 > 1/6) then
        ion_vx_mm = ion_vx_mm - photon_mom / _amu_mass_per_charge
      elseif (random_number_emit_1 <= 3/6 and random_number_emit_1 > 2/6) then
        ion_vy_mm = ion_vy_mm + photon_mom / _amu_mass_per_charge
      elseif (random_number_emit_1 <= 4/6 and random_number_emit_1 > 3/6) then
        ion_vy_mm = ion_vy_mm - photon_mom / _amu_mass_per_charge
      elseif (random_number_emit_1 <= 5/6 and random_number_emit_1 > 4/6) then
        ion_vz_mm = ion_vz_mm + photon_mom / _amu_mass_per_charge
      elseif (random_number_emit_1 <= 1 and random_number_emit_1 > 5/6) then
        ion_vz_mm = ion_vz_mm - photon_mom / _amu_mass_per_charge
      end


    end

    local random_number_2 = rand()

    if (R_z_per_dt > random_number_2) then

      ion_vz_mm = ion_vz_mm + sign(z_0,z,ion_vz_mm) * photon_mom / _amu_mass_per_charge

      local KE_z = speed_to_ke(ion_vz_mm,100)
      -- file_z:write(ion_time_of_flight, ',') 
      -- file_z:write(KE_z, ',') 
      -- file_z:write(R_z_per_dt, '\n')

      local random_number_emit_2 = rand()

      if (random_number_emit_2 <= 1/6 and random_number_emit_2 > 0) then
        ion_vx_mm = ion_vx_mm + photon_mom / _amu_mass_per_charge
      elseif (random_number_emit_2 <= 2/6 and random_number_emit_2 > 1/6) then
        ion_vx_mm = ion_vx_mm - photon_mom / _amu_mass_per_charge
      elseif (random_number_emit_2 <= 3/6 and random_number_emit_2 > 2/6) then
        ion_vy_mm = ion_vy_mm + photon_mom / _amu_mass_per_charge
      elseif (random_number_emit_2 <= 4/6 and random_number_emit_2 > 3/6) then
        ion_vy_mm = ion_vy_mm - photon_mom / _amu_mass_per_charge
      elseif (random_number_emit_2 <= 5/6 and random_number_emit_2 > 4/6) then
        ion_vz_mm = ion_vz_mm + photon_mom / _amu_mass_per_charge
      elseif (random_number_emit_2 <= 1 and random_number_emit_2 > 5/6) then
        ion_vz_mm = ion_vz_mm - photon_mom / _amu_mass_per_charge
      end


    end

    


  end


end

--------------------------------------------------------------------------------

-- SIMION segment called by SIMION to override time-step size on each time-step.
function segment.tstep_adjust()

   -- Keep time step size below some fraction of the RF period.
   -- See "Time Step Size" comments.
   --ion_time_step = min(ion_time_step, 0.1*1E+6/frequency_hz)  -- X usec
   --ion_time_step = 0.01*1E+6/frequency_hz
     ion_time_step = 0.00739 -- mus

end

--------------------------------------------------------------------------------

function segment.terminate_run()
  --print('num hits on test plan:', num_hits)
  file:close()
  file_x:close()
  file_y:close()
  file_z:close()

    -- calculate/display emittance
    --print("Num particles = " .. #y)
    --local x_emit, norm_x_emit = compute_x_emittance(x, xprime, vx, vz)
    --local y_emit, norm_y_emit = compute_y_emittance(y, yprime, vy, vz)
    --local total_emit = x_emit * y_emit
    --local total_norm_emit = norm_x_emit * norm_y_emit
    --print("Beam Emittance = " .. total_emit .. " mm * mrad (Normalized = " .. total_norm_emit .. ")")

  --------------------------------------------------------------------------------
    --local KE = compute_KE(vx_v2, vy_v2, vz_v2)
    --print(KE)
end

function segment.terminate()

  --file:close()
  --file_v2:close()

end
