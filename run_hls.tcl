############################################################
## Vitis HLS TCL Script for Heston Monte Carlo
## Target: Zynq Ultrascale+ ZCU102 (xczu9eg-ffvb1156-2-e)
############################################################

# Create new project
open_project heston_mc_project

# Set top function
set_top heston_mc_kernel

# Add design files
add_files heston_mc.cpp -cflags "-std=c++11"
add_files heston_mc.h -cflags "-std=c++11"
add_files Mersenne_Twister.cpp -cflags "-std=c++11"
add_files Mersenne_Twister.h -cflags "-std=c++11"

# Add testbench
add_files -tb heston_mc_tb.cpp -cflags "-std=c++11"

# Create solution targeting ZCU102
open_solution "solution1" -flow_target vivado

# Set ZCU102 part
set_part {xczu9eg-ffvb1156-2-e}

# Create clock with 5ns period (200MHz)
create_clock -period 5 -name default

# Set synthesis directives
config_compile -name_max_length 80 -pipeline_loops 64

# Interface synthesis options
config_interface -m_axi_latency 64 -m_axi_max_widen_bitwidth 512

# RTL configuration
config_rtl -reset all -reset_async -reset_level low

# Configure scheduling
config_schedule -effort high -enable_dsp_full_reg

# Run C simulation
csim_design -clean

# Run synthesis
csynth_design

# Run C/RTL co-simulation (optional, takes longer)
# Uncomment the following line to run co-simulation
# cosim_design -trace_level all

# Export IP for Vivado integration
export_design -format ip_catalog -description "Heston Monte Carlo Option Pricing" -vendor "user" -version "1.0"

# Exit
exit