library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package regionizer_data is 
    type particle is record
        pt : signed(15 downto 0);
        eta : signed(9 downto 0);
        phi : signed(9 downto 0);
        rest : std_logic_vector(27 downto 0);
    end record;
    type glbparticle is record
        pt : signed(15 downto 0);
        eta : signed(11 downto 0);
        phi : signed(10 downto 0);
        rest : std_logic_vector(24 downto 0);
    end record;



    subtype word64 is std_logic_vector(63 downto 0);
    subtype word65 is std_logic_vector(64 downto 0);

    function particle_to_w64(p : particle) return word64;
    function w64_to_particle(d : word64) return particle;
    function null_particle return particle;
    function glbparticle_to_w64(p : glbparticle) return word64;
    function w64_to_glbparticle(d : word64) return glbparticle;
    function null_glbparticle return glbparticle;

    type particles is array(natural range <>) of particle;
    type w64s      is array(natural range <>) of word64;
    type w65s      is array(natural range <>) of word65;
    type glbparticles is array(natural range <>) of glbparticle;

    constant PHI_SHIFT_INT : natural := 160; -- 2*pi/9, size of a phi nonant, track finder sector or fiducial part of one PF region
    constant PHI_SHIFT : signed(9 downto 0) := to_signed(160, 10); -- 2*pi/9, size of a phi nonant, track finder sector or fiducial part of one PF region
    constant PHI_BORDER : signed(9 downto 0) := to_signed(58, 10); -- 0.25 (0.30 would be 69) 
    constant PHI_MARGIN_POS : signed(9 downto 0) := to_signed(+(160/2-58), 10);  -- half-width of fiducial MINUS border (half-size of gap are between sector N and sector N+2)
    constant PHI_MARGIN_NEG : signed(9 downto 0) := to_signed(-(160/2-58), 10);  -- same but with negative sign
    constant PHI_HALFWIDTH_POS : signed(9 downto 0) := to_signed(+(160/2+58), 10); -- half size of a full region (fiducial PLUS border)
    constant PHI_HALFWIDTH_NEG : signed(9 downto 0) := to_signed(-(160/2+58), 10);  
    constant ETA_HALFWIDTH_POS : signed(9 downto 0) := to_signed(+(230/2+58), 10); -- half size of a full region (fiducial PLUS border)
    constant ETA_HALFWIDTH_NEG : signed(9 downto 0) := to_signed(-(230/2+58), 10);  

    constant PHI_CALOSHIFT    : signed(9 downto 0) := to_signed( 480,        10);  -- 2*pi/3, size of an HGCal sector
    constant PHI_CALOSHIFT1   : signed(9 downto 0) := to_signed( 320,        10);  -- 2*pi/3 - 2*pi/9, distance between center of hgcal sector 1 and pf region 1 = 2 * size of a phi nonant
    constant PHI_CALOEDGE_POS : signed(9 downto 0) := to_signed(+(480/2-58), 10);  -- +(half-size of calo sector)-border
    constant PHI_CALOEDGE_NEG : signed(9 downto 0) := to_signed(-(480/2-58), 10);  -- -(half-size of calo sector)+border

    constant PHI_MPI : signed(11 downto 0) := to_signed(-PHI_SHIFT_INT*9/2, 12);  -- same but with negative sign
    constant PHI_2PI : signed(11 downto 0) := to_signed( PHI_SHIFT_INT*9,   12);  -- same but with negative sign

    constant PFII : natural := 6;
    constant PFII240 : natural := 4;
    constant NPFREGIONS : natural := 9;

    constant NTKSECTORS : natural := 9;
    constant NTKFIBERS : natural := 2;
    constant NTKFIFOS : natural := NTKFIBERS*3;
    constant NTKSORTED : natural := 30;
    constant NTKSTREAM : natural := (NTKSORTED+PFII240-1)/PFII240;

    constant NCALOSECTORS : natural := 3;
    constant NCALOFIBERS : natural := 4;
    constant NCALOFIFO0 : natural := NCALOFIBERS;
    constant NCALOFIFO12 : natural := 2*NCALOFIBERS;
    constant NCALOFIFOS : natural := NCALOFIFO0+2*NCALOFIFO12;
    constant NCALOSORTED : natural := 20;
    constant NCALOSTREAM : natural := (NCALOSORTED+PFII240-1)/PFII240;

    constant NMUFIBERS : natural := 2;
    constant NMUSORTED : natural := 4;
    constant NMUSTREAM : natural := (NMUSORTED+PFII240-1)/PFII240;

    constant NPFTOT :     natural := NTKSORTED + NCALOSORTED + NMUSORTED;
    constant NPFSTREAM :  natural := (NPFTOT+PFII240-1)/PFII240;

    constant NPUPPI :     natural := NTKSORTED + NCALOSORTED;
    constant NPUPPIFINALSORTED : natural := 18;

    constant TKDELAY : natural := 2;
    constant MUDELAY   : natural := 3;
    constant CALODELAY : natural := 1;

    constant TDEMUX_FACTOR      : natural := 3;
    constant TDEMUX_NTKFIBERS   : natural := NTKFIBERS*3/2/TDEMUX_FACTOR;
    constant TDEMUX_NCALOFIBERS : natural := NCALOFIBERS;
    constant TDEMUX_NMUFIBERS   : natural := 1;

end package;

package body regionizer_data is
    function particle_to_w64(p : particle) return word64 is
        variable ret : word64;
    begin
        ret( 9 downto  0) := std_logic_vector(p.eta);
        ret(19 downto 10) := std_logic_vector(p.phi);
        ret(47 downto 32) := std_logic_vector(p.pt);
        ret(31 downto 20) := (p.rest(11 downto  0));
        ret(63 downto 48) := (p.rest(27 downto 12));
        return ret;
    end particle_to_w64;

    function w64_to_particle(d : word64) return particle is
        variable ret : particle;
    begin
        ret.eta := signed(d( 9 downto  0));
        ret.phi := signed(d(19 downto 10));
        ret.pt  := signed(d(47 downto 32));
        ret.rest(11 downto  0) := (d(31 downto 20));
        ret.rest(27 downto 12) := (d(63 downto 48));
        return ret;
    end w64_to_particle;

    function null_particle return particle is
        variable ret : particle;
    begin
        ret.eta := to_signed(0, ret.eta'length);
        ret.phi := to_signed(0, ret.phi'length);
        ret.pt  := to_signed(0,  ret.pt'length);
        ret.rest := (others => '0');
        return ret;
    end null_particle;

    function glbparticle_to_w64(p : glbparticle) return word64 is
        variable ret : word64;
    begin
        ret(11 downto  0) := std_logic_vector(p.eta);
        ret(22 downto 12) := std_logic_vector(p.phi);
        ret(47 downto 32) := std_logic_vector(p.pt);
        ret(31 downto 23) := (p.rest( 8 downto 0));
        ret(63 downto 48) := (p.rest(24 downto 9));
        return ret;
    end glbparticle_to_w64;

    function w64_to_glbparticle(d : word64) return glbparticle is
        variable ret : glbparticle;
    begin
        ret.eta := signed(d(11 downto  0));
        ret.phi := signed(d(22 downto 12));
        ret.pt  := signed(d(47 downto 32));
        ret.rest( 8 downto 0) := (d(31 downto 23));
        ret.rest(24 downto 9) := (d(63 downto 48));
        return ret;
    end w64_to_glbparticle;

    function null_glbparticle return glbparticle is
        variable ret : glbparticle;
    begin
        ret.eta := to_signed(0, ret.eta'length);
        ret.phi := to_signed(0, ret.phi'length);
        ret.pt  := to_signed(0,  ret.pt'length);
        ret.rest := (others => '0');
        return ret;
    end null_glbparticle;

end regionizer_data;


