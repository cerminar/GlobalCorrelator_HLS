-- ==============================================================
-- RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and SystemC
-- Version: 2018.1
-- Copyright (C) 1986-2018 Xilinx, Inc. All Rights Reserved.
-- 
-- ===========================================================

library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity dr2_int_cap_14_s is
port (
    ap_ready : OUT STD_LOGIC;
    eta1_V : IN STD_LOGIC_VECTOR (9 downto 0);
    phi1_V : IN STD_LOGIC_VECTOR (9 downto 0);
    eta2_V : IN STD_LOGIC_VECTOR (9 downto 0);
    phi2_V : IN STD_LOGIC_VECTOR (9 downto 0);
    ap_return : OUT STD_LOGIC_VECTOR (12 downto 0) );
end;


architecture behav of dr2_int_cap_14_s is 
    constant ap_const_logic_1 : STD_LOGIC := '1';
    constant ap_const_boolean_1 : BOOLEAN := true;
    constant ap_const_lv11_0 : STD_LOGIC_VECTOR (10 downto 0) := "00000000000";
    constant ap_const_lv32_7 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000000111";
    constant ap_const_lv32_A : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000001010";
    constant ap_const_lv4_0 : STD_LOGIC_VECTOR (3 downto 0) := "0000";
    constant ap_const_lv15_1063 : STD_LOGIC_VECTOR (14 downto 0) := "001000001100011";
    constant ap_const_lv13_1063 : STD_LOGIC_VECTOR (12 downto 0) := "1000001100011";
    constant ap_const_logic_0 : STD_LOGIC := '0';

    signal lhs_V_fu_48_p1 : STD_LOGIC_VECTOR (10 downto 0);
    signal rhs_V_fu_52_p1 : STD_LOGIC_VECTOR (10 downto 0);
    signal r_V_fu_56_p2 : STD_LOGIC_VECTOR (10 downto 0);
    signal tmp_fu_62_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_s_fu_68_p2 : STD_LOGIC_VECTOR (10 downto 0);
    signal deta_V_fu_74_p3 : STD_LOGIC_VECTOR (10 downto 0);
    signal lhs_V_1_fu_86_p1 : STD_LOGIC_VECTOR (10 downto 0);
    signal rhs_V_1_fu_90_p1 : STD_LOGIC_VECTOR (10 downto 0);
    signal r_V_1_fu_94_p2 : STD_LOGIC_VECTOR (10 downto 0);
    signal tmp_5_fu_100_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_6_fu_106_p2 : STD_LOGIC_VECTOR (10 downto 0);
    signal dphi_V_fu_112_p3 : STD_LOGIC_VECTOR (10 downto 0);
    signal tmp_7_fu_124_p2 : STD_LOGIC_VECTOR (10 downto 0);
    signal tmp_1248_fu_130_p4 : STD_LOGIC_VECTOR (3 downto 0);
    signal deta2_V_fu_184_p2 : STD_LOGIC_VECTOR (13 downto 0);
    signal dphi2_V_fu_191_p2 : STD_LOGIC_VECTOR (13 downto 0);
    signal rhs_V_2_fu_149_p1 : STD_LOGIC_VECTOR (14 downto 0);
    signal lhs_V_2_fu_146_p1 : STD_LOGIC_VECTOR (14 downto 0);
    signal icmp_fu_140_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal r_V_2_fu_152_p2 : STD_LOGIC_VECTOR (14 downto 0);
    signal val_assign_fu_158_p3 : STD_LOGIC_VECTOR (14 downto 0);
    signal tmp_9_fu_166_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_1249_fu_172_p1 : STD_LOGIC_VECTOR (12 downto 0);
    signal deta2_V_fu_184_p0 : STD_LOGIC_VECTOR (10 downto 0);
    signal deta_V_cast2_fu_82_p1 : STD_LOGIC_VECTOR (13 downto 0);
    signal deta2_V_fu_184_p1 : STD_LOGIC_VECTOR (10 downto 0);
    signal dphi2_V_fu_191_p0 : STD_LOGIC_VECTOR (10 downto 0);
    signal dphi_V_cast1_fu_120_p1 : STD_LOGIC_VECTOR (13 downto 0);
    signal dphi2_V_fu_191_p1 : STD_LOGIC_VECTOR (10 downto 0);

    component mp7wrapped_pfalgoeOg IS
    generic (
        ID : INTEGER;
        NUM_STAGE : INTEGER;
        din0_WIDTH : INTEGER;
        din1_WIDTH : INTEGER;
        dout_WIDTH : INTEGER );
    port (
        din0 : IN STD_LOGIC_VECTOR (10 downto 0);
        din1 : IN STD_LOGIC_VECTOR (10 downto 0);
        dout : OUT STD_LOGIC_VECTOR (13 downto 0) );
    end component;



begin
    mp7wrapped_pfalgoeOg_U272 : component mp7wrapped_pfalgoeOg
    generic map (
        ID => 1,
        NUM_STAGE => 1,
        din0_WIDTH => 11,
        din1_WIDTH => 11,
        dout_WIDTH => 14)
    port map (
        din0 => deta2_V_fu_184_p0,
        din1 => deta2_V_fu_184_p1,
        dout => deta2_V_fu_184_p2);

    mp7wrapped_pfalgoeOg_U273 : component mp7wrapped_pfalgoeOg
    generic map (
        ID => 1,
        NUM_STAGE => 1,
        din0_WIDTH => 11,
        din1_WIDTH => 11,
        dout_WIDTH => 14)
    port map (
        din0 => dphi2_V_fu_191_p0,
        din1 => dphi2_V_fu_191_p1,
        dout => dphi2_V_fu_191_p2);




    ap_ready <= ap_const_logic_1;
    ap_return <= 
        tmp_1249_fu_172_p1 when (tmp_9_fu_166_p2(0) = '1') else 
        ap_const_lv13_1063;
    deta2_V_fu_184_p0 <= deta_V_cast2_fu_82_p1(11 - 1 downto 0);
    deta2_V_fu_184_p1 <= deta_V_cast2_fu_82_p1(11 - 1 downto 0);
        deta_V_cast2_fu_82_p1 <= std_logic_vector(IEEE.numeric_std.resize(signed(deta_V_fu_74_p3),14));

    deta_V_fu_74_p3 <= 
        r_V_fu_56_p2 when (tmp_fu_62_p2(0) = '1') else 
        tmp_s_fu_68_p2;
    dphi2_V_fu_191_p0 <= dphi_V_cast1_fu_120_p1(11 - 1 downto 0);
    dphi2_V_fu_191_p1 <= dphi_V_cast1_fu_120_p1(11 - 1 downto 0);
        dphi_V_cast1_fu_120_p1 <= std_logic_vector(IEEE.numeric_std.resize(signed(dphi_V_fu_112_p3),14));

    dphi_V_fu_112_p3 <= 
        r_V_1_fu_94_p2 when (tmp_5_fu_100_p2(0) = '1') else 
        tmp_6_fu_106_p2;
    icmp_fu_140_p2 <= "1" when (tmp_1248_fu_130_p4 = ap_const_lv4_0) else "0";
        lhs_V_1_fu_86_p1 <= std_logic_vector(IEEE.numeric_std.resize(signed(phi1_V),11));

    lhs_V_2_fu_146_p1 <= std_logic_vector(IEEE.numeric_std.resize(unsigned(deta2_V_fu_184_p2),15));
        lhs_V_fu_48_p1 <= std_logic_vector(IEEE.numeric_std.resize(signed(eta1_V),11));

    r_V_1_fu_94_p2 <= std_logic_vector(signed(lhs_V_1_fu_86_p1) - signed(rhs_V_1_fu_90_p1));
    r_V_2_fu_152_p2 <= std_logic_vector(unsigned(rhs_V_2_fu_149_p1) + unsigned(lhs_V_2_fu_146_p1));
    r_V_fu_56_p2 <= std_logic_vector(signed(lhs_V_fu_48_p1) - signed(rhs_V_fu_52_p1));
        rhs_V_1_fu_90_p1 <= std_logic_vector(IEEE.numeric_std.resize(signed(phi2_V),11));

    rhs_V_2_fu_149_p1 <= std_logic_vector(IEEE.numeric_std.resize(unsigned(dphi2_V_fu_191_p2),15));
        rhs_V_fu_52_p1 <= std_logic_vector(IEEE.numeric_std.resize(signed(eta2_V),11));

    tmp_1248_fu_130_p4 <= tmp_7_fu_124_p2(10 downto 7);
    tmp_1249_fu_172_p1 <= val_assign_fu_158_p3(13 - 1 downto 0);
    tmp_5_fu_100_p2 <= "1" when (signed(r_V_1_fu_94_p2) > signed(ap_const_lv11_0)) else "0";
    tmp_6_fu_106_p2 <= std_logic_vector(unsigned(ap_const_lv11_0) - unsigned(r_V_1_fu_94_p2));
    tmp_7_fu_124_p2 <= (dphi_V_fu_112_p3 or deta_V_fu_74_p3);
    tmp_9_fu_166_p2 <= "1" when (unsigned(val_assign_fu_158_p3) < unsigned(ap_const_lv15_1063)) else "0";
    tmp_fu_62_p2 <= "1" when (signed(r_V_fu_56_p2) > signed(ap_const_lv11_0)) else "0";
    tmp_s_fu_68_p2 <= std_logic_vector(unsigned(ap_const_lv11_0) - unsigned(r_V_fu_56_p2));
    val_assign_fu_158_p3 <= 
        r_V_2_fu_152_p2 when (icmp_fu_140_p2(0) = '1') else 
        ap_const_lv15_1063;
end behav;
