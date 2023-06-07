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
adjustable frequency_hz             = 1.1E6 --1.1E6  -- RF frequency of quad (Hz)

adjustable q_u = 0.620
adjustable a_u = 0.122

--------------------------------------------------------------------------------

adjustable T1_min = 17
adjustable T1_max = 22
adjustable T2_min = 19
adjustable T2_max = 23

adjustable V_endcap_min = 70
adjustable V_endcap_max = 130

adjustable V_add = 4

local i = 0
local j = 0
local k = 0

adjustable i_clone = 0
adjustable j_clone = 0
adjustable k_clone = 0

--------------------------------------------------------------------------------

-- various arrays to store variables on each particle
local x = {}        -- x positions (mm)
local xprime = {}   -- x' (radians)
local y = {}        -- y positions (mm)
local yprime = {}   -- y' (radians)
local vx = {}       -- x-velocity (mm/usec)
local vy = {}       -- y-velocity (mm/usec)
local vz = {}       -- z-velocity (mm/usec)

--------------------------------------------------------------------------------

function segment.flym() -- Called at the beginning of every flym
  sim_trajectory_image_control = 1 -- Don't preserve trajectories

for k = V_endcap_min, V_endcap_max, 10 do
    for i = T1_min, T1_max, 0.1 do
      for j = T2_min, T2_max, 0.1 do
        i_clone = i
        j_clone = j
        k_clone = k
        print('T1 = ' .. i_clone .. ' mus and T2 = ' .. j_clone .. ' mus')
        run()
      end
    end
  end
end


--------------------------------------------------------------------------------

function compute_x_emittance(x, xprime, vx, vz)
  -- Compute average of all numbers in given array.
  -- Returns 0 if array contains zero elements.
  function average(array)
      local result = 0
      for _,a in ipairs(array) do result = result + a end
      if #array ~= 0 then result = result / #array end
      return result
  end

  -- Compute various averages for emittance.
  local x_ave = average(x)
  local xprime_ave = average(xprime)
  local t = {}; for n = 1,#x do t[n] = (x[n] - x_ave)^2 end
  local dx2_ave = average(t)
  local t = {}; for n = 1,#x do t[n] = (xprime[n] - xprime_ave)^2 end
  local dxprime2_ave = average(t)
  local t = {}; for n = 1,#x do t[n] = (x[n]-x_ave)*(xprime[n]-xprime_ave) end
  local dx_dxprime_ave = average(t)

  -- Compute emittance from averages, in correct units.
  local m = dx2_ave * dxprime2_ave - dx_dxprime_ave^2
  if m < 0 then m = 0 end      -- safety on numerical roundoff
  local x_emit = sqrt(m) * 1000  -- (mm * mrad)

  -- Compute average speed for normalized emittance.
  local vx_avg = average(vx)
  local vz_avg = average(vz)
  local v_avg = sqrt(vx_avg^2 + vz_avg^2)
  --FIX: or this:
  --local t = {}; for n = 1,#x do t[n] = sqrt(vx[n]^2 + vx[n]^2) end
  --local v_avg = average(t)

  -- compute normalized emittance from averages
  local c = 300000                    -- speed of light (mm/usec)
  local beta = v_avg / c              -- relativistic beta
  local gamma = 1 / sqrt(1 - beta^2)  -- relativistic gamma
  local norm_x_emit = beta * gamma * x_emit

  return x_emit, norm_x_emit
end



function compute_y_emittance(y, yprime, vy, vz)
    -- Compute average of all numbers in given array.
    -- Returns 0 if array contains zero elements.
    function average(array)
        local result = 0
        for _,a in ipairs(array) do result = result + a end
        if #array ~= 0 then result = result / #array end
        return result
    end

    -- Compute various averages for emittance.
    local y_ave = average(y)
    local yprime_ave = average(yprime)
    local t = {}; for n = 1,#y do t[n] = (y[n] - y_ave)^2 end
    local dy2_ave = average(t)
    local t = {}; for n = 1,#y do t[n] = (yprime[n] - yprime_ave)^2 end
    local dyprime2_ave = average(t)
    local t = {}; for n = 1,#y do t[n] = (y[n]-y_ave)*(yprime[n]-yprime_ave) end
    local dy_dyprime_ave = average(t)

    -- Compute emittance from averages, in correct units.
    local m = dy2_ave * dyprime2_ave - dy_dyprime_ave^2
    if m < 0 then m = 0 end      -- safety on numerical roundoff
    local y_emit = sqrt(m) * 1000  -- (mm * mrad)

    -- Compute average speed for normalized emittance.
    local vy_avg = average(vy)
    local vz_avg = average(vz)
    local v_avg = sqrt(vy_avg^2 + vz_avg^2)
    --FIX: or this:
    --local t = {}; for n = 1,#y do t[n] = sqrt(vx[n]^2 + vy[n]^2) end
    --local v_avg = average(t)

    -- compute normalized emittance from averages
    local c = 300000                    -- speed of light (mm/usec)
    local beta = v_avg / c              -- relativistic beta
    local gamma = 1 / sqrt(1 - beta^2)  -- relativistic gamma
    local norm_y_emit = beta * gamma * y_emit

    return y_emit, norm_y_emit
end

--------------------------------------------------------------------------------

-- Note: Using circular rods, the radius of the rods themselves
-- should optimally be approximately 1.1487 * r_0.

-- Temporary variables used internally.
local scaled_rf  -- a factor used in the RF component
local omega      -- frequency_hz (reexpressed in units of radians/usec)
local theta      -- phase_angle_deg (reexpressed in units of radians)
local last_pe_update = 0.0 -- last potential energy surface update time (usec)

-- SIMION segment called by SIMION to set adjustable electrode voltages
-- in the current potential array instance.
-- NOTE: this is called frequently, multiple times per time-step (by
-- Runge-Kutta), so performance concerns here can be important.

--------------------------------------------------------------------------------

local num_hits = 0

local TP = simion.import 'testplanelib.lua'

local test = TP(0,0,0.11,0,0,-1,
  -- example of function to call on reaching test plane.
  function()
    --mark()
    --print('In test plane: n = ' .. ion_number .. ' z = ' .. ion_pz_mm)
    -- ion_splat = 1  -- optionally splat particle in test plane
    num_hits = num_hits + 1  -- optionally count hits on test plane

    -- store variables for emittance calculation
    local particle_count = #y + 1
    x[particle_count] = ion_px_mm    -- store x position (mm)
    y[particle_count] = ion_py_mm    -- store y position (mm)
    vx[particle_count] = ion_vx_mm   -- store x-velocity (mm/usec)
    vy[particle_count] = ion_vy_mm   -- store y-velocity (mm/usec)
    vz[particle_count] = ion_vz_mm   -- store y-velocity (mm/usec)
    xprime[particle_count] = ion_vx_mm / ion_vz_mm  -- store ~tan(theta) (rad)
    yprime[particle_count] = ion_vy_mm / ion_vz_mm  -- store ~tan(theta) (rad)

    --x[particle_count] = ion_px_mm    -- store x position (mm)
    --xprime[particle_count] = ion_vx_mm / ion_vz_mm  -- store ~tan(theta) (rad)
    -- FIX? or this: yprime[particle_count] = atan2(ion_vy_mm, ion_vx_mm)

  end

)

--------------------------------------------------------------------------------



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

  adj_elect01 = tempvolts - dcvolts + dcvolts *0
  adj_elect02 = - dcvolts *0

  adj_elect03 = tempvolts - dcvolts + dcvolts *0
  adj_elect04 = - dcvolts *0

  adj_elect05 = tempvolts - dcvolts + dcvolts *0
  adj_elect06 = - dcvolts *0

  adj_elect07 = tempvolts - dcvolts + dcvolts *0
  adj_elect08 = - dcvolts *0

  adj_elect09 = tempvolts - dcvolts + dcvolts *0
  adj_elect10 = - dcvolts *0

  adj_elect11 = tempvolts - dcvolts + dcvolts *0
  adj_elect12 = - dcvolts *0

  adj_elect13 = tempvolts - dcvolts + dcvolts *0
  adj_elect14 = - dcvolts *0

  if (ion_time_of_flight <= i_clone) then
    adj_elect15 = tempvolts - dcvolts + dcvolts
    adj_elect16 = - dcvolts
  elseif (ion_time_of_flight > i_clone and ion_time_of_flight <= j_clone) then
    adj_elect15 = tempvolts - dcvolts + dcvolts
    adj_elect16 = - dcvolts
  elseif (ion_time_of_flight > j_clone) then
    adj_elect15 = tempvolts - dcvolts + dcvolts + V_add
    adj_elect16 = - dcvolts + V_add
  end

  adj_elect17 = tempvolts - dcvolts + dcvolts
  adj_elect18 = - dcvolts

  if (ion_time_of_flight <= i_clone) then
    adj_elect19 = tempvolts - dcvolts + dcvolts
    adj_elect20 = - dcvolts
  elseif (ion_time_of_flight > i_clone and ion_time_of_flight <= j_clone) then
    adj_elect19 = tempvolts - dcvolts + dcvolts + k_clone/2
    adj_elect20 = - dcvolts + k_clone/2 -- Ruben's suggestion
  elseif (ion_time_of_flight > j_clone) then
    adj_elect19 = tempvolts - dcvolts + dcvolts + V_add
    adj_elect20 = - dcvolts + V_add
  end


end




--------------------------------------------------------------------------------

local printcounter = 0

-- SIMION segment called by SIMION after every time-step.
function segment.other_actions()
  -- Update potential energy surface display periodically.
  -- The performance overhead of this in non-PE views is only a few percent.
  -- NOTE: the value inside abs(...) can be negative when a new ion is flown.
  if abs(ion_time_of_flight - last_pe_update) >= pe_update_each_usec then
    last_pe_update = ion_time_of_flight
    sim_update_pe_surface = 1    -- Request a PE surface display update.
  end


  printcounter = printcounter + 1 -- To record data properly
  if (printcounter == 20) then
    --mark()
    printcounter = 0
  end

  if (ion_splat ~= 0 ) then -- Unstable trajectory
    --print("Ion : " .. ion_number .. " , UNSTABLE")
    print("UNSTABLE")
    mark()
  elseif (ion_time_of_flight > 200) then
    --print("Ion : " .. ion_number .. " , STABLE")
    print("STABLE")
    mark()
    ion_splat = - 4 -- Ion killed
  end

  test.other_actions()

  --if (abs(ion_vz_mm) <= 0.05) then
    --print(ion_time_of_flight)
  --end

end

--------------------------------------------------------------------------------

-- SIMION segment called by SIMION to override time-step size on each time-step.
function segment.tstep_adjust()
   -- Keep time step size below some fraction of the RF period.
   -- See "Time Step Size" comments.
   ion_time_step = min(ion_time_step, 0.1*1E+6/frequency_hz)  -- X usec

   test.tstep_adjust()
   --- saus
end

--------------------------------------------------------------------------------



  --function segment.tstep_adjust()
    --test.tstep_adjust()
  --end

  --function segment.other_actions()
    --test.other_actions()
  --end

function segment.terminate_run()
  --print('num hits on test plan:', num_hits)

  -- calculate/display emittance
  --print("Num particles = " .. #y)
  local x_emit, norm_x_emit = compute_x_emittance(x, xprime, vx, vz)
  local y_emit, norm_y_emit = compute_y_emittance(y, yprime, vy, vz)
  local total_emit = x_emit * y_emit
  local total_norm_emit = norm_x_emit * norm_y_emit
  --print("Beam Emittance = " .. total_emit .. " mm * mrad (Normalized = " .. total_norm_emit .. ")")

end
