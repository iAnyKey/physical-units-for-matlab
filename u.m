classdef u < handle
% u  Physical units.
%
%   If the Physical Units Toolbox is on your MATLAB path, there is nothing to
%   initialize, add to your workspace, or pass to functions - simply
%   multiply/divide by u.(unitName) to attach physical units to a variable. For
%   example, to define a speed using a supported unit: carSpeed = 100 * u.kph.
%   Or, define a speed with an unsupported unit as a combination of supported
%   units: snailSpeed = 20 * u.m/u.week.
%
%   Calling u by itself will display all available units in u.
%
%   Variables with physical units attached are of the class DimVar
%   ("dimenensioned variable"). Math operations performed on dimensioned
%   variables will automatically perform dimensional analysis and can create new
%   units or cancel units and return a normal variable.
%
%   When displaying variables with units or using them in plot, etc., the units
%   used for display will be, in order if available and valid: 1) per-variable
%   preferred display units, 2) units listed in displayUnits, 3) a combination
%   of fundamental base units (mass, length, time, temperature, ...). To set (or
%   clear) custom display units for a given variable, see the scd function. To
%   customize the displayUnits list, see displayUnits. For more advanced
%   customization of the base units themselves, see baseUnitSystem.
%
%   Display customization is set by calls to displayUnits and/or baseUnitSystem
%   (either function files or variables in the base workspace). Tailor
%   preferences for a specific project by defining these variables at the top of
%   a script (before any units are called) or placing unique versions of the
%   files in a project's directory. Be sure to clear the class when changing
%   projects or else the old customizations will remain in effect, e.g.
%       <a href="matlab:clear u">clear u</a> or 
%       <a href="matlab:clear classes">clear classes</a>.
%
%   Some MATLAB functions won't accept variables with physical units (DimVars).
%   Most of the time displayingvalue, which returns value in terms of preferred
%   display units, will be the appropriate tool, but there is also double, which
%   returns the value in terms of base units, and u2num.
%
%   Example 1: Shaft power.
%       rotationSpeed = 2500 * u.rpm;
%       torque = 95 * str2u('ft-lbf');  % Use alternate string-based definition.
%       power = rotationSpeed * torque; % Returns variable with units of power.
%       horsePower = power / u.hp;      % Convert/cancel units.
%
%   Example 2: Unit conversion.
%       100 * u.acre/u.ha;  % Convert 100 acres to hectares.
%       u.st/u.kg;          % Return conversion factor for stone to kilos.
% 
%   Example 3: Custom display.
%       fieldSize = 3*u.sqkm
%       % Muliplies and divides remove any per-variable custom display units.
%       % It's nice to display in the units that make sense for that variable.
%       rate = 3.7*u.acre/u.day
%       rate = scd(rate,'sqm/hr')
%       timeNeeded = fieldSize/rate
%       timeNeeded = scd(timeNeeded,'month')
% 
%   See also displayUnits, baseUnitSystem, scd, clear, displayingvalue,
%   DimVar.double, u2num, str2u, symunit,
%     dispdisp - http://www.mathworks.com/matlabcentral/fileexchange/48637.

%   Copyright Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715
%   github.com/sky-s/physical-units-for-matlab

properties (Hidden, Constant = true)
    %% User-defined base and display:
    % Establishes base unit system and preferences based on baseUnitSystem and
    % displayUnits.
     baseUnitSystem =    evalEverythere_('baseUnitSystem')
     dispUnits =         evalEverythere_('displayUnits')
     coreUnits =         buildCoreUnits(evalEverythere_('u.baseUnitSystem'))
end
properties (Constant = true)

    %% Core units:
    baseNames = evalEverythere_("u.baseUnitSystem(:,1)")';
                             

    m               = evalEverythere_("scd(u.coreUnits.m,'m')") % meter
    kg              = evalEverythere_("scd(u.coreUnits.kg,'kg')") % kilogram
    s               = evalEverythere_("scd(u.coreUnits.s,'s')") % second
    A               = evalEverythere_("scd(u.coreUnits.A,'A')") % ampere
    K               = evalEverythere_("scd(u.coreUnits.K,'K')") % kelvin (°C = °K-273.15)
    mol             = evalEverythere_("scd(u.coreUnits.mol,'mol')") % mole
    cd              = evalEverythere_("scd(u.coreUnits.cd,'cd')") % candela
    bit             = evalEverythere_("scd(u.coreUnits.bit,'bit')") % bit
    currency        = evalEverythere_("scd(u.coreUnits.currency,'currency')") % currency
    unit            = evalEverythere_("scd(u.coreUnits.unit,'unit')") % user unit
    
    %% SI defining constants
    D_Cs            = evalEverythere_("scd(9192631770/u.s,'D_Cs')") % hyperfine transition frequency of Cs
    c               = evalEverythere_("scd(299792458*u.m/u.s,'c')") % speed of light in vacuum
    c_0             = evalEverythere_("scd(u.c,'c_0')") % speed of light in vacuum
    lightSpeed      = evalEverythere_("scd(u.c,'lightSpeed')") 
    speedOfLight    = evalEverythere_("scd(u.c,'speedOfLight')") 
    h_c             = evalEverythere_("scd(6.62607015e-34*u.kg*u.m^2/u.s,'h_c')") % Planck constant
    PlanckConstant  = evalEverythere_("scd(u.h_c,'PlanckConstant')") 
    e               = evalEverythere_("scd(1.602176634e-19*u.A*u.s,'e')") % elementary charge
    elementaryCharge= evalEverythere_("scd(u.e,'elementaryCharge')") 
    k               = evalEverythere_("scd(1.380649e-23*u.kg*u.m^2/u.s^2/u.K,'k')") % Boltzmann constant
    k_B             = evalEverythere_("scd(u.k,'k_B')") % Boltzmann constant
    BoltzmannConstant = evalEverythere_("scd(u.k,'BoltzmannConstant')")
    N_A             = evalEverythere_("scd(6.02214076e23/u.mol,'N_A')") % Avogadro constant
    NA              = evalEverythere_("scd(u.N_A,'NA')") % Avogadro constant
    AvogadroConstant= evalEverythere_("scd(u.N_A,'AvogadroConstant')")
    K_cd            = evalEverythere_("scd(683*u.cd/(u.kg*u.m^2/u.s^3),'K_cd')") % luminous efficacy of 540 THz radiation
    
    h_bar           = evalEverythere_("scd(u.h_c/(2*pi),'h_bar')") % Dirac constant
    DiracConstant   = evalEverythere_("scd(u.h_bar,'DiracConstant')")
    
    %% SI prefixes
    yotta = 1e24; % SI prefix yotta, Y 
    zetta = 1e21; % SI prefix zetta, Z
    exa = 1e18; % SI prefix exa, E
    peta = 1e15; % SI prefix peta, P
    tera = 1e12; % SI prefix tera, T
    giga = 1e9; % SI prefix giga, G
    mega = 1e6; % SI prefix mega, M
    kilo = 1e3; % SI prefix kilo, k
    hecto = 1e2; % SI prefix hecto, h
    deka = 1e1; % SI prefix deka, da
    deci = 1e-1; % SI prefix deci, d
    centi = 1e-2; % SI prefix centi, c
    milli = 1e-3; % SI prefix milli, m
    micro = 1e-6; % SI prefix micro, µ
    nano = 1e-9; % SI prefix nano, n
    pico = 1e-12; % SI prefix pico, p
    femto = 1e-15; % SI prefix femto, f
    atto = 1e-18; % SI prefix atto, a
    zepto = 1e-21; % SI prefix zepto, z
    yocto = 1e-24; % SI prefix yocto, y
    
    %% Important constants with uncertainty
    alpha                   = 7.2973525693e-3 % fine-structure constant
    fine_structureConstant  = evalEverythere_("u.alpha")
    SommerfeldConstant      = evalEverythere_("u.alpha")
    G                       = evalEverythere_("scd(6.67430e-11*u.m^3/u.kg/u.s^2,'G')") % gravitational constant
    gravitationalConstant   = evalEverythere_("scd(u.G,'gravitationalConstant')") 
    m_e                     = evalEverythere_("scd(9.1093837015e-31*u.kg,'m_e')") % electron rest mass
    electronMass            = evalEverythere_("scd(u.m_e,'electronMass')")
    m_p                     = evalEverythere_("scd(1.67262192369e-27*u.kg,'m_e')") % proton rest mass

    %% Derived units list:
    % References:
    % http://physics.nist.gov/cuu/Constants/index.html
    % http://www.translatorscafe.com/unit-converter
    % http://en.wikipedia.org
    % http://www.efunda.com/units/index.cfm
    
    %---- length ----

    meter               = evalEverythere_("scd(u.m,'meter')") 
    km                  = evalEverythere_("scd(1e3*u.m,'km')") % kilometer
    kilometer           = evalEverythere_("scd(u.km,'kilometer')") 
    dm                  = evalEverythere_("scd(1e-1*u.m,'dm')") % decimeter
    decimeter           = evalEverythere_("scd(u.dm,'decimeter')") 
    cm                  = evalEverythere_("scd(1e-2*u.m,'cm')") % centimeter
    centimeter          = evalEverythere_("scd(u.cm,'centimeter')") 
    mm                  = evalEverythere_("scd(1e-3*u.m,'mm')") % millimeter
    millimeter          = evalEverythere_("scd(u.mm,'millimeter')") 
    um                  = evalEverythere_("scd(1e-6*u.m,'um')") % micrometer
    micrometer          = evalEverythere_("scd(u.um,'micrometer')") 
    micron              = evalEverythere_("scd(u.um,'micron')") % micron
    nm                  = evalEverythere_("scd(1e-9*u.m,'nm')") % nanometer
    nanometer           = evalEverythere_("scd(u.nm,'nanometer')") 
    pm                  = evalEverythere_("scd(1e-12*u.m,'pm')") % picometer
    picometer           = evalEverythere_("scd(u.pm,'picometer')") 
    fm                  = evalEverythere_("scd(1e-15*u.m,'fm')") % femtometer
    femtometer          = evalEverythere_("scd(u.fm,'femtometer')") 
    fermi               = evalEverythere_("scd(u.fm,'fermi')") % fermi
    Ao                  = evalEverythere_("scd(1e-10*u.m,'Ao')") % ångström
    ang                 = evalEverythere_("scd(u.Ao,'ang')") % ångström
    angstrom            = evalEverythere_("scd(u.ang,'angstrom')") 
    angstroem           = evalEverythere_("scd(u.ang,'angstroem')") 
    a0                  = evalEverythere_("scd(u.h_bar/(u.m_e*u.c*u.alpha),'a0')") % Bohr radius
    a_0                 = evalEverythere_("scd(u.a0,'a_0')") % Bohr radius
    BohrRadius          = evalEverythere_("scd(u.a0,'BohrRadius')")
    lP                  = evalEverythere_("scd(u.h_bar/(u.mP*u.c),'lP')") % Planck length
    PlanckLength        = evalEverythere_("scd(u.lP,'PlanckLength')") 
    xu                  = evalEverythere_("scd(1.0021e-13*u.m,'xu')") % x unit
    xUnit               = evalEverythere_("scd(u.xu,'xUnit')") 
    xu_Cu               = evalEverythere_("scd(1.00207697e-13*u.m,'xu_Cu')") % x unit (copper)
    xUnit_copper        = evalEverythere_("scd(u.xu_Cu,'xUnit_copper')") 
    xu_Mo               = evalEverythere_("scd(1.00209952e-13*u.m,'xu_Mo')") % x unit (molybdenum)
    xUnit_molybdenum    = evalEverythere_("scd(u.xu_Mo,'xUnit_molybdenum')") 
    in                  = evalEverythere_("scd(2.54*u.cm,'in')") % inch
    inch                = evalEverythere_("scd(u.in,'inch')") 
    mil                 = evalEverythere_("scd(1e-3*u.in,'mil')") % mil
    line                = evalEverythere_("scd(u.in/10,'line')") % line
    pica                = evalEverythere_("scd(u.in/6,'pica')") % pica
    point               = evalEverythere_("scd(u.pica/12,'point')") % point
    hand                = evalEverythere_("scd(4*u.in,'hand')") % hand
    span                = evalEverythere_("scd(9*u.in,'span')") % span
    smoot               = evalEverythere_("scd(67*u.in,'smoot')") % smoot
    ft                  = evalEverythere_("scd(12*u.in,'ft')") % foot
    foot                = evalEverythere_("scd(u.ft,'foot')") 
    ft_US               = evalEverythere_("scd(1200/3937*u.m,'ft_US')") % US survey foot
    foot_US             = evalEverythere_("scd(u.ft_US,'foot_US')") % US survey foot
    kft                 = evalEverythere_("scd(1e3*u.ft,'kft')") % kilofoot
    kilofoot            = evalEverythere_("scd(u.kft,'kilofoot')") 
    FL                  = evalEverythere_("scd(100*u.ft,'FL')") % flight level
    flightLevel         = evalEverythere_("scd(u.FL,'flightLevel')") 
    yd                  = evalEverythere_("scd(3*u.ft,'yd')") % yard
    yard                = evalEverythere_("scd(u.yd,'yard')") 
    ftm                 = evalEverythere_("scd(6*u.ft,'ftm')") % fathom
    fathom              = evalEverythere_("scd(u.ftm,'fathom')")
    li                  = evalEverythere_("scd(0.66*u.ft,'li')") % link
    link                = evalEverythere_("scd(u.li,'link')") 
    rod                 = evalEverythere_("scd(5.5*u.yd,'rod')") % rod
    ch                  = evalEverythere_("scd(66*u.ft,'ch')") % chain
    chain               = evalEverythere_("scd(u.ch,'chain')") 
    fur                 = evalEverythere_("scd(220*u.yd,'fur')") % furlong
    furlong             = evalEverythere_("scd(u.fur,'furlong')") 
    mi                  = evalEverythere_("scd(5280*u.ft,'mi')") % mile
    mile                = evalEverythere_("scd(u.mi,'mile')") 
    mi_US               = evalEverythere_("scd(6336/3937*u.km,'mi_US')") % US survey mile
    mile_US             = evalEverythere_("scd(u.mi_US,'mile_US')") % US survey mile
    nmi                 = evalEverythere_("scd(1852*u.m,'nmi')") % nautical mile
    NM                  = evalEverythere_("scd(u.nmi,'NM')") % nautical mile
    inm                 = evalEverythere_("scd(u.nmi,'inm')") % nautical mile
    nauticalMile        = evalEverythere_("scd(u.nmi,'nauticalMile')") 
    nm_UK               = evalEverythere_("scd(6080*u.ft,'nm_UK')") % Imperial nautical mile
    nmile               = evalEverythere_("scd(u.nm_UK,'nmile')") % Imperial nautical mile
    dataMile            = evalEverythere_("scd(6000*u.ft,'dataMile')") 
    au                  = evalEverythere_("scd(149597870.7*u.km,'au')") % astronomical unit
    astronomicalUnit    = evalEverythere_("scd(u.au,'astronomicalUnit')") 
    pc                  = evalEverythere_("scd(648000/pi*u.au,'pc')") % parsec
    parsec              = evalEverythere_("scd(u.pc,'parsec')") 

    %---- reciprocal length ----

    dpt                 = evalEverythere_("scd(1/u.m,'dpt')") % diopter
    diopter             = evalEverythere_("scd(u.dpt,'diopter')") 
    R_inf               = evalEverythere_("scd(u.alpha/(4*pi*u.a0),'R_inf')") % Rydberg constant
    RydbergConstant     = evalEverythere_("scd(u.R_inf,'RydbergConstant')") 

    %---- area ----

    ft2                 = evalEverythere_("scd(u.ft^2,'ft²')") % square foot
    sqft                = evalEverythere_("scd(u.ft2,'sqft')") % square foot
    square              = evalEverythere_("scd(100*u.sqft,'square')") % square
    ha                  = evalEverythere_("scd(10000*u.m^2,'ha')") % hectare
    hectare             = evalEverythere_("scd(u.ha,'hectare')") 
    a                   = evalEverythere_("scd(100*u.m^2,'a')") % are
    are                 = evalEverythere_("scd(u.a,'are')") 
    ac                  = evalEverythere_("scd(43560*u.sqft,'ac')") % acre
    acre                = evalEverythere_("scd(u.ac,'acre')") 
    ro                  = evalEverythere_("scd(1/4*u.acre,'ro')") % rood
    rood                = evalEverythere_("scd(u.ro,'rood')") 
    twp                 = evalEverythere_("scd(36*u.mi^2,'twp')") % township
    township            = evalEverythere_("scd(u.twp,'township')") 
    circ_mil            = evalEverythere_("scd(pi/4*u.mil^2,'circ_mil')") % circular mil
    circularMil         = evalEverythere_("scd(u.circ_mil,'circularMil')") 
    circ_inch           = evalEverythere_("scd(pi/4*u.in^2,'circ_inch')") % circular inch
    circularInch        = evalEverythere_("scd(u.circ_inch,'circularInch')") 
    b                   = evalEverythere_("scd(100*u.fm^2,'b')") % barn
    barn                = evalEverythere_("scd(u.b,'barn')") 
    sqin                = evalEverythere_("scd(u.in^2,'sqin')") % square inch
    squareInch          = evalEverythere_("scd(u.sqin,'squareInch')")
    sqmil               = evalEverythere_("scd(u.mil^2,'sqmil')") % square mil
    squareMil           = evalEverythere_("scd(u.sqmil,'squareMil')")
    sqmi                = evalEverythere_("scd(u.mi^2,'sqmi')") % square mile
    squareMile          = evalEverythere_("scd(u.sqmi,'squareMile')")
    sqnmi               = evalEverythere_("scd(u.nmi^2,'sqnmi')") % square nautical mile
    squareNauticalMile  = evalEverythere_("scd(u.sqnmi,'squareNauticalMile')")
    m2                  = evalEverythere_("scd(u.m^2,'m²')") % square meter
    sqm                 = evalEverythere_("scd(u.m^2,'sqm')") % square meter
    squareMeter         = evalEverythere_("scd(u.sqm,'squareMeter')")
    sqkm                = evalEverythere_("scd(u.km^2,'sqkm')") % square kilometer
    squareKilometer     = evalEverythere_("scd(u.sqkm,'squareKilometer')")
    sqcm                = evalEverythere_("scd(u.cm^2,'sqcm')") % square centimeter
    squareCentimeter    = evalEverythere_("scd(u.sqcm,'squareCentimeter')")
    sqmm                = evalEverythere_("scd(u.mm^2,'sqmm')") % square millimeter
    squareMillimeter    = evalEverythere_("scd(u.sqmm,'squareMillimeter')")
    sqdm                = evalEverythere_("scd(u.dm^2,'sqdm')") % square decimeter
    squareDecimeter     = evalEverythere_("scd(u.sqdm,'squareDecimeter')")

    %---- volume ----

    m3                  = evalEverythere_("scd(u.m^3,'m³')") % cubic meter
    cc                  = evalEverythere_("scd(u.cm^3,'cc')") % cubic centimeter
    cubicCentimeter     = evalEverythere_("scd(u.cc,'cubicCentimeter')") 
    L                   = evalEverythere_("scd(1000*u.cc,'L')") % liter
    l                   = evalEverythere_("scd(u.L,'l')") % liter
    liter               = evalEverythere_("scd(u.L,'liter')") 
    dL                  = evalEverythere_("scd(100*u.cc,'dL')") % deciliter
    dl                  = evalEverythere_("scd(u.dL,'dl')") % deciliter
    deciliter           = evalEverythere_("scd(u.dl,'deciliter')") 
    cL                  = evalEverythere_("scd(10*u.cc,'cL')") % centiliter
    cl                  = evalEverythere_("scd(u.cL,'cl')") % centiliter
    centiliter          = evalEverythere_("scd(u.cl,'centiliter')") 
    mL                  = evalEverythere_("scd(u.cc,'mL')") % milliliter
    ml                  = evalEverythere_("scd(u.mL,'ml')") % milliliter
    milliliter          = evalEverythere_("scd(u.ml,'milliliter')") 
    uL                  = evalEverythere_("scd(u.mm^3,'uL')") % microliter
    ul                  = evalEverythere_("scd(u.uL,'ul')") % microliter
    microliter          = evalEverythere_("scd(u.ul,'microliter')") 
    kL                  = evalEverythere_("scd(u.m^3,'kL')") % kiloliter
    kl                  = evalEverythere_("scd(u.kL,'kl')") % kiloliter
    kiloliter           = evalEverythere_("scd(u.kl,'kiloliter')") 
    cuin                = evalEverythere_("scd(16.387064*u.mL,'cuin')") % cubic inch
    cubicInch           = evalEverythere_("scd(u.cuin,'cubicInch')")
    ft3                 = evalEverythere_("scd(u.ft^3,'ft³')") % cubic foot
    FBM                 = evalEverythere_("scd(u.sqft*u.in,'FBM')") % board foot
    boardFoot           = evalEverythere_("scd(u.FBM,'boardFoot')") 
    CCF                 = evalEverythere_("scd(100*u.ft3,'CCF')") % centum cubic feet
    HCF                 = evalEverythere_("scd(u.CCF,'HCF')") % hundred cubic feet
    MCF                 = evalEverythere_("scd(1000*u.ft3,'MCF')") % mille cubic feet
    MMCF                = evalEverythere_("scd(1e6*u.ft3,'MMCF')") % mille mille (1e6) cubic feet
    BCF                 = evalEverythere_("scd(1e9*u.ft3,'BCF')") % billion cubic feet
    TMC                 = evalEverythere_("scd(u.BCF,'TMC')") % thousand million cubic feet
    TCF                 = evalEverythere_("scd(1e12*u.ft3,'TCF')") % trillion cubic feet
    gal                 = evalEverythere_("scd(231*u.cuin,'gal')") % gallon (US)
    gallon              = evalEverythere_("scd(u.gal,'gallon')") 
    gal_UK              = evalEverythere_("scd(4.54609*u.l,'gal_UK')") % UK imperial gallon
    igal                = evalEverythere_("scd(u.gal_UK,'igal')") % UK imperial gallon
    quart               = evalEverythere_("scd(u.gal/4,'quart')") % US quart
    qt_UK               = evalEverythere_("scd(u.gal_UK/4,'qt_UK')") % British imperial quart
    liq_qt              = evalEverythere_("scd(u.quart,'liq_qt')") % US quart
    pint                = evalEverythere_("scd(u.quart/2,'pint')") % US pint
    pint_UK             = evalEverythere_("scd(u.qt_UK/2,'pint_UK')") % British imperial pint
    liq_pt              = evalEverythere_("scd(u.pint,'liq_pt')") % US pint
    cup                 = evalEverythere_("scd(u.pint/2,'cup')") % US cup
    floz                = evalEverythere_("scd(u.cup/8,'floz')") % US fluid ounce
    fluidOunce          = evalEverythere_("scd(u.floz,'fluidOunce')") % US fluid ounce
    floz_UK             = evalEverythere_("scd(u.gal_UK/160,'floz_UK')") % British imperial fluid ounce
    Tbls                = evalEverythere_("scd(u.floz/2,'Tbls')") % US tablespoon
    tablespoon          = evalEverythere_("scd(u.Tbls,'tablespoon')") % US tablespoon
    tsp                 = evalEverythere_("scd(u.Tbls/3,'tsp')") % US teaspoon
    teaspoon            = evalEverythere_("scd(u.tsp,'teaspoon')") % US teaspoon
    acft                = evalEverythere_("scd(u.acre*u.ft,'acft')") % acre-foot
    acre_foot           = evalEverythere_("scd(u.acft,'acre_foot')") 
    acin                = evalEverythere_("scd(u.acre*u.in,'acin')") % acre-inch
    acre_inch           = evalEverythere_("scd(u.acin,'acre_inch')") 
    bbl                 = evalEverythere_("scd(7056*u.in^3,'bbl')") % US customary dry barrel
    barrel              = evalEverythere_("scd(u.bbl,'barrel')") 
    fldr                = evalEverythere_("scd(u.floz/8,'fldr')") % US customary fluid dram
    fluidDram           = evalEverythere_("scd(u.fldr,'fluidDram')") 
    fldr_UK             = evalEverythere_("scd(u.floz_UK/8,'fldr_UK')") % British imperial fluid drachm (dram)
    minim               = evalEverythere_("scd(u.fldr/60,'minim')") % US customary minim
    minim_UK            = evalEverythere_("scd(u.fldr_UK/60,'minim_UK')") % British imperial minim
    gill                = evalEverythere_("scd(4*u.floz,'gill')") % US customary fluid gill
    gill_UK             = evalEverythere_("scd(u.gal_UK/32,'gill_UK')") % British imperial gill

    %---- acceleration ----

    g0                  = evalEverythere_("scd(9.80665*u.m/u.s^2,'g0')") % standard gravity
    g_0                 = evalEverythere_("scd(u.g0,'g_0')") % standard gravity
    gn                  = evalEverythere_("scd(u.g0,'gn')") % standard gravity
    g_n                 = evalEverythere_("scd(u.g0,'g_n')") % standard gravity
    gee                 = evalEverythere_("scd(u.g0,'gee')") % standard gravity
    standardGravity     = evalEverythere_("scd(u.g0,'standardGravity')") 
    Gal                 = evalEverythere_("scd(u.cm/u.s^2,'Gal')") % gal

    %---- force ----

    N                   = evalEverythere_("scd(u.kg*u.m/u.s^2,'N')") % newton
    newton              = evalEverythere_("scd(u.N,'newton')") 
    kN                  = evalEverythere_("scd(1e3*u.N,'kN')") % kilonewton
    kilonewton          = evalEverythere_("scd(u.kN,'kilonewton')") 
    MN                  = evalEverythere_("scd(1e6*u.N,'MN')") % meganewton
    meganewton          = evalEverythere_("scd(u.MN,'meganewton')") 
    mN                  = evalEverythere_("scd(1e-3*u.N,'mN')") % millinewton
    millinewton         = evalEverythere_("scd(u.mN,'millinewton')") 
    uN                  = evalEverythere_("scd(1e-6*u.N,'uN')") % micronewton
    micronewton         = evalEverythere_("scd(u.uN,'micronewton')") 
    dyn                 = evalEverythere_("scd(1e-5*u.N,'dyn')") % dyne
    dyne                = evalEverythere_("scd(u.dyn,'dyne')") 
    lbf                 = evalEverythere_("scd(4.4482216152605*u.N,'lbf')") % pound force
    lb_f                = evalEverythere_("scd(u.lbf,'lb_f')") % pound force
    poundForce          = evalEverythere_("scd(u.lbf,'poundForce')") 
    kip                 = evalEverythere_("scd(1000*u.lbf,'kip')") % kip
    kilopoundForce      = evalEverythere_("scd(u.kip,'kilopoundForce')") 
    kgf                 = evalEverythere_("scd(u.kg*u.g0,'kgf')") % kilogram force
    kg_f                = evalEverythere_("scd(u.kgf,'kg_f')") % kilogram force
    kilogramForce       = evalEverythere_("scd(u.kgf,'kilogramForce')") 
    kp                  = evalEverythere_("scd(u.kgf,'kp')") % kilopond
    kilopond            = evalEverythere_("scd(u.kp,'kilopond')") 
    p                   = evalEverythere_("scd(u.kp/1000,'p')") % pond
    pond                = evalEverythere_("scd(u.p,'pond')") 
    sn                  = evalEverythere_("scd(u.kN,'sn')") % sthène
    sthene              = evalEverythere_("scd(u.sn,'sthene')") 

    %---- mass ----

    kilogram            = evalEverythere_("scd(u.kg,'kilogram')")
    g                   = evalEverythere_("scd(1e-3*u.kg,'g')") % gram
    gram                = evalEverythere_("scd(u.g,'gram')") 
    mg                  = evalEverythere_("scd(1e-3*u.gram,'mg')") % milligram
    milligram           = evalEverythere_("scd(u.mg,'milligram')") 
    ug                  = evalEverythere_("scd(1e-6*u.gram,'ug')") % microgram
    microgram           = evalEverythere_("scd(u.ug,'microgram')") 
    Mg                  = evalEverythere_("scd(1e6*u.gram,'Mg')") % Megagram/metric tonne
    Megagram            = evalEverythere_("scd(u.Mg,'Megagram')") 
    Gg                  = evalEverythere_("scd(1e9*u.gram,'Gg')") % Gigagram
    Gigagram            = evalEverythere_("scd(u.Gg,'Gigagram')") 
    t                   = evalEverythere_("scd(1000*u.kg,'t')") % metric tonne
    tonne               = evalEverythere_("scd(u.t,'tonne')") % metric ton
    Mt                  = evalEverythere_("scd(1e6*u.t,'Mt')") % metric megatonne
    megatonne           = evalEverythere_("scd(u.Mt,'megatonne')") 
    lbm                 = evalEverythere_("scd(0.45359237*u.kg,'lbm')") % pound mass
    lb_m                = evalEverythere_("scd(u.lbm,'lb_m')") % pound mass
    poundMass           = evalEverythere_("scd(u.lbm,'poundMass')") 
    lb                  = evalEverythere_("scd(u.lbm,'lb')") % pound mass
    pound               = evalEverythere_("scd(u.lb,'pound')") 
    tn                  = evalEverythere_("scd(2000*u.lbm,'tn')") % US customary short ton
    ton                 = evalEverythere_("scd(u.tn,'ton')") % US customary short ton
    ton_UK              = evalEverythere_("scd(2240*u.lbm,'ton_UK')") % British imperial ton
    st                  = evalEverythere_("scd(14*u.lbm,'st')") % stone
    stone               = evalEverythere_("scd(u.st,'stone')") 
    cwt                 = evalEverythere_("scd(100*u.lbm,'cwt')") % US customary short hundredweight
    hundredweight       = evalEverythere_("scd(u.cwt,'hundredweight')") 
    cwt_UK              = evalEverythere_("scd(8*u.stone,'cwt_UK')") % British imperial short hundredweight
    quarter             = evalEverythere_("scd(u.cwt_UK/4,'quarter')") % British imperial quarter
    slug                = evalEverythere_("scd(u.lbf/(u.ft/u.s^2),'slug')") % slug
    slinch              = evalEverythere_("scd(u.lbf/(u.in/u.s^2),'slinch')") 
    blob                = evalEverythere_("scd(u.slinch,'blob')")
    oz                  = evalEverythere_("scd(u.lbm/16,'oz')") % ounce
    ounce               = evalEverythere_("scd(u.oz,'ounce')") 
    dr                  = evalEverythere_("scd(u.oz/16,'dr')") % dram
    dram                = evalEverythere_("scd(u.dr,'dram')") 
    gr                  = evalEverythere_("scd(u.lbm/7000,'gr')") % grain
    grain               = evalEverythere_("scd(u.gr,'grain')") 
    ct                  = evalEverythere_("scd(200*u.mg,'ct')") % carat
    carat               = evalEverythere_("scd(u.ct,'carat')") 
    amu                 = evalEverythere_("scd(1.66053906660e-27*u.kg,'amu')") % atomic mass unit
    atomicMassUnit      = evalEverythere_("scd(u.amu,'atomicMassUnit')") 
    Da                  = evalEverythere_("scd(u.amu,'Da')") % atomic mass unit
    dalton              = evalEverythere_("scd(u.Da,'dalton')") 
    mu                  = evalEverythere_("scd(u.amu,'mu')") % atomic mass unit
    mP                  = evalEverythere_("scd(sqrt(u.h_bar*u.c/u.G),'mP')") % Planck mass
    PlanckMass          = evalEverythere_("scd(u.mP,'PlanckMass')") 
    mug                 = evalEverythere_("scd(u.kgf/(u.m/u.s^2),'mug')") % metric slug
    metricSlug          = evalEverythere_("scd(u.mug,'metricSlug')") 
    hyl                 = evalEverythere_("scd(u.mug,'hyl')") % hyl
    par                 = evalEverythere_("scd(u.mug,'par')") % par
    TMU                 = evalEverythere_("scd(u.mug,'TMU')") % technische Masseneinheit
    technischeMasseneinheit = evalEverythere_("scd(u.TMU,'technischeMasseneinheit')") 
    glug                = evalEverythere_("scd(u.g*u.g0/(u.cm/u.s^2),'glug')") 

    %---- more force ----

    pdl                 = evalEverythere_("scd(u.lbm*u.ft/u.s^2,'pdl')") % poundal
    poundal             = evalEverythere_("scd(u.pdl,'poundal')") 
    gf                  = evalEverythere_("scd(u.gram*u.g0,'gf')") % gram force
    g_f                 = evalEverythere_("scd(u.gf,'g_f')") %gram force
    gramForce           = evalEverythere_("scd(u.gf,'gramForce')") 
    ozf                 = evalEverythere_("scd(u.oz*u.g0,'ozf')") % ounce force
    oz_f                = evalEverythere_("scd(u.ozf,'oz_f')") % ounce force
    ounceForce          = evalEverythere_("scd(u.ozf,'ounceForce')") 
    tonf                = evalEverythere_("scd(u.tn*u.g0,'tonf')") % short ton force
    ton_f               = evalEverythere_("scd(u.tonf,'ton_f')") % short ton force
    tonForce            = evalEverythere_("scd(u.tonf,'tonForce')")

    %---- mass per length ----

    den                 = evalEverythere_("scd(u.gram/(9*u.km),'den')") % denier
    denier              = evalEverythere_("scd(u.den,'denier')") 
    tex                 = evalEverythere_("scd(u.gram/u.km,'tex')") % tex
    dtex                = evalEverythere_("scd(u.tex/10,'dtex')") % decitex
    decitex             = evalEverythere_("scd(u.dtex,'decitex')")

    %---- time ----

    second              = evalEverythere_("scd(u.s,'second')") 
    sec                 = evalEverythere_("scd(u.s,'sec')") % second
    ms                  = evalEverythere_("scd(1e-3*u.s,'ms')") % millisecond
    millisecond         = evalEverythere_("scd(u.ms,'millisecond')") 
    us                  = evalEverythere_("scd(1e-6*u.s,'us')") % microsecond
    microsecond         = evalEverythere_("scd(u.us,'microsecond')") 
    ns                  = evalEverythere_("scd(1e-9*u.s,'ns')") % nanosecond
    nanosecond          = evalEverythere_("scd(u.ns,'nanosecond')") 
    ps                  = evalEverythere_("scd(1e-12*u.s,'ps')") % picosecond
    picosecond          = evalEverythere_("scd(u.ps,'picosecond')") 
    fs                  = evalEverythere_("scd(1e-15*u.s,'fs')") % femtosecond
    femtosecond         = evalEverythere_("scd(u.fs,'femtosecond')") 
    tP                  = evalEverythere_("scd(u.lP/u.c,'tP')") % Planck time
    PlanckTime          = evalEverythere_("scd(u.tP,'PlanckTime')") 
    min                 = evalEverythere_("scd(60*u.s,'min')") % minute
    minute              = evalEverythere_("scd(u.min,'minute')") 
    h                   = evalEverythere_("scd(60*u.min,'h')") % hour
    hr                  = evalEverythere_("scd(u.h,'hr')") % hour
    hrs                 = evalEverythere_("scd(u.h,'hrs')")
    hour                = evalEverythere_("scd(u.hr,'hour')") 
    d                   = evalEverythere_("scd(24*u.hr,'d')") % day
    day                 = evalEverythere_("scd(u.d,'day')") % day
    days                = evalEverythere_("scd(u.d,'days')")
    week                = evalEverythere_("scd(7*u.day,'week')") % week
    fortnight           = evalEverythere_("scd(2*u.week,'fortnight')") % fortnight
    month_30            = evalEverythere_("scd(30*u.day,'month_30')") % 30-day month
    yr                  = evalEverythere_("scd(365.25*u.day,'yr')") % julian year
    y                   = evalEverythere_("scd(u.yr,'y')") % julian year
    year                = evalEverythere_("scd(u.yr,'year')") % julian year
    yrs                 = evalEverythere_("scd(u.yr,'yrs')") % julian year
    year_julian         = evalEverythere_("scd(u.year,'year_julian')") % julian year
    year_360            = evalEverythere_("scd(360*u.day,'year_360')") % 360-day year
    year_Tropical       = evalEverythere_("scd(365.24219*u.day,'year_Tropical')") % tropical year
    year_Gregorian      = evalEverythere_("scd(365.2425*u.day,'year_Gregorian')") % gregorian year
    month               = evalEverythere_("scd(u.yr/12,'month')") % 1/12th julian year
    flick               = evalEverythere_("scd(u.s/705600000,'flick')") 

    %---- frequency ----

    Hz                  = evalEverythere_("scd(1/u.s,'Hz')") % hertz (NB: incompatible with angle and angular velocity)
    hertz               = evalEverythere_("scd(u.Hz,'hertz')") 
    kHz                 = evalEverythere_("scd(1e3*u.Hz,'kHz')") % kilohertz
    kilohertz           = evalEverythere_("scd(u.kHz,'kilohertz')") 
    MHz                 = evalEverythere_("scd(1e6*u.Hz,'MHz')") % megahertz
    megahertz           = evalEverythere_("scd(u.MHz,'megahertz')") 
    GHz                 = evalEverythere_("scd(1e9*u.Hz,'GHz')") % gigahertz
    gigahertz           = evalEverythere_("scd(u.GHz,'gigahertz')") 
    THz                 = evalEverythere_("scd(1e12*u.Hz,'THz')") % terahertz
    terahertz           = evalEverythere_("scd(u.THz,'terahertz')") 
    Bd                  = evalEverythere_("scd(1/u.s,'Bd')") % baud
    baud                = evalEverythere_("scd(u.Bd,'baud')") 

    %---- energy ----

    Nm                  = evalEverythere_("scd(u.N*u.m,'Nm')") % newton-meter
    newton_meter        = evalEverythere_("scd(u.Nm,'newton_meter')") 
    J                   = evalEverythere_("scd(u.Nm,'J')") % joule
    joule               = evalEverythere_("scd(u.J,'joule')") 
    kJ                  = evalEverythere_("scd(1e3*u.J,'kJ')") % kilojoule
    kilojoule           = evalEverythere_("scd(u.kJ,'kilojoule')") 
    MJ                  = evalEverythere_("scd(1e6*u.J,'MJ')") % megajoule
    megajoule           = evalEverythere_("scd(u.MJ,'megajoule')") 
    GJ                  = evalEverythere_("scd(1e9*u.J,'GJ')") % gigajoule
    gigajoule           = evalEverythere_("scd(u.GJ,'gigajoule')") 
    TJ                  = evalEverythere_("scd(1e12*u.J,'TJ')") % terajoule
    terajoule           = evalEverythere_("scd(u.TJ,'terajoule')")
    mJ                  = evalEverythere_("scd(1e-3*u.J,'mJ')") % millijoule
    millijoule          = evalEverythere_("scd(u.mJ,'millijoule')") 
    uJ                  = evalEverythere_("scd(1e-6*u.J,'uJ')") % microjoule
    microjoule          = evalEverythere_("scd(u.uJ,'microjoule')") 
    nJ                  = evalEverythere_("scd(1e-9*u.J,'nJ')") % nanojoule
    nanojoule           = evalEverythere_("scd(u.nJ,'nanojoule')") 
    eV                  = evalEverythere_("scd(u.e/u.C*u.J,'eV')") % electronvolt
    electronvolt        = evalEverythere_("scd(u.eV,'electronvolt')") 
    BTU                 = evalEverythere_("scd(1055.06*u.J,'BTU')") % British thermal unit (ISO)
    Btu                 = evalEverythere_("scd(u.BTU,'Btu')") % British thermal unit (ISO)
    kBtu                = evalEverythere_("scd(1e3*u.Btu,'kBtu')") % kiloBtu
    MMBtu               = evalEverythere_("scd(1e6*u.Btu,'MMBtu')") % million Btu
    britishThermalUnit  = evalEverythere_("scd(u.Btu,'britishThermalUnit')") 
    Btu_IT              = evalEverythere_("scd(1055.05585*u.J,'Btu_IT')") % British thermal unit (International Table)
    Btu_th              = evalEverythere_("scd(1054.3503*u.J,'Btu_th')") % British thermal unit (thermochemical)
    kpm                 = evalEverythere_("scd(u.kp*u.m,'kpm')") % kilopond-meter
    kilopond_meter      = evalEverythere_("scd(u.kpm,'kilopond_meter')") 
    Ws                  = evalEverythere_("scd(u.J,'Ws')") % watt-second
    watt_second         = evalEverythere_("scd(u.Ws,'watt_second')") 
    kWh                 = evalEverythere_("scd(3.6e6*u.J,'kWh')") % kilowatt-hour
    kilowatt_hour       = evalEverythere_("scd(u.kWh,'kilowatt_hour')") 
    Wh                  = evalEverythere_("scd(3.6e3*u.J,'Wh')") % watt-hour
    watt_hour           = evalEverythere_("scd(u.Wh,'watt_hour')") 
    cal                 = evalEverythere_("scd(4.1868*u.J,'cal')") % calorie (International Table)
    calorie             = evalEverythere_("scd(u.cal,'calorie')") 
    cal_IT              = evalEverythere_("scd(u.cal,'cal_IT')") % calorie (International Table)
    cal_4               = evalEverythere_("scd(4.204*u.J,'cal_4')") % calorie (4°C)
    cal_15              = evalEverythere_("scd(4.1855*u.J,'cal_15')") % calorie (15°C)
    cal_20              = evalEverythere_("scd(4.182*u.J,'cal_20')") % calorie (20°C)
    cal_mean            = evalEverythere_("scd(4.190*u.J,'cal_mean')") % calorie (mean)
    cal_th              = evalEverythere_("scd(4.184*u.J,'cal_th')") % calorie (thermochemical)
    tonTnt              = evalEverythere_("scd(u.cal_th*1e9,'tonTnt')") % ton of TNT
    kcal                = evalEverythere_("scd(1e3*u.cal,'kcal')") % kilocalorie
    kilocalorie         = evalEverythere_("scd(u.kcal,'kilocalorie')") 
    kcal_IT             = evalEverythere_("scd(1e3*u.cal_IT,'kcal_IT')") % kilocalorie (International Table)
    Cal                 = evalEverythere_("scd(u.kcal,'Cal')") % large calorie / food calorie
    foodCalorie         = evalEverythere_("scd(u.Cal,'foodCalorie')") 
    largeCalorie        = evalEverythere_("scd(u.Cal,'largeCalorie')") 
    kcal_4              = evalEverythere_("scd(1e3*u.cal_4,'kcal_4')") % kilocalorie (4°C)
    kcal_15         	= evalEverythere_("scd(1e3*u.cal_15,'kcal_15')") % kilocalorie (15°C)
    kcal_20             = evalEverythere_("scd(1e3*u.cal_20,'kcal_20')") % kilocalorie (20°C)
    kcal_mean           = evalEverythere_("scd(1e3*u.cal_mean,'kcal_mean')") % kilocalorie (mean)
    kcal_th             = evalEverythere_("scd(1e3*u.cal_th,'kcal_th')") % kilocalorie (thermochemical)
    erg                 = evalEverythere_("scd(1e-7*u.J,'erg')") % en.wikipedia.org/wiki/Erg
    E_h                 = evalEverythere_("scd(u.alpha^2*u.m_e*u.c^2,'E_h')") % Hartree energy
    Ha                  = evalEverythere_("scd(u.E_h,'Ha')") % hartree
    hartree             = evalEverythere_("scd(u.Ha,'hartree')") 
    thm                 = evalEverythere_("scd(1e5*u.BTU,'thm')") % therm
    therm               = evalEverythere_("scd(u.thm,'therm')") 
    quad                = evalEverythere_("scd(1e15*u.BTU,'quad')") % quad

    %---- temperature ----
    % For reference: °C = °K-273.15; °F = °R-459.67.

    kelvin              = evalEverythere_("scd(u.K,'kelvin')") 
    degK                = evalEverythere_("scd(u.K,'°K')") % degrees kelvin
    R                   = evalEverythere_("scd(u.K*5/9,'R')") % rankine (°F = °R-459.67)
    rankine             = evalEverythere_("scd(u.R,'rankine')") 
    degR                = evalEverythere_("scd(u.R,'°R')") % degrees rankine
    degC                = evalEverythere_("scd(OffsetDimVar(u.K,273.15*u.K),'°C')") % Celcius
    Celcius             = evalEverythere_("scd(u.degC,'Celcius')");
    centigrade          = evalEverythere_("scd(u.degC,'centigrade')");
    degF                = evalEverythere_("scd(OffsetDimVar(u.R,459.67*u.R),'°F')") % Fahrenheit
    Fahrenheit          = evalEverythere_("scd(u.degF,'Fahrenheit')");
    % Réaumur
    % Rømer
    DeltaK              = evalEverythere_("scd(u.K,'DeltaK')") % kelvin (relative temperature)
    DeltadegC           = evalEverythere_("scd(u.K,'Delta°C')") % celsius (relative, °C = °K-273.15)
    DeltadegR           = evalEverythere_("scd(u.R,'Delta°R')") % rankine (relative temperature)
    DeltadegF           = evalEverythere_("scd(u.R,'Delta°F')") % fahrenheit (relative, °F = °R-459.67)
    mK                  = evalEverythere_("scd(1e-3*u.K,'mK')") % millikelvin
    millikelvin         = evalEverythere_("scd(u.mK,'millikelvin')") 
    uK                  = evalEverythere_("scd(1e-6*u.K,'uK')") % microkelvin
    microkelvin         = evalEverythere_("scd(u.uK,'microkelvin')") 
    nK                  = evalEverythere_("scd(1e-9*u.K,'nK')") % nanokelvin
    nanokelvin          = evalEverythere_("scd(u.nK,'nanokelvin')") 
    TP                  = evalEverythere_("scd(u.mP*u.c^2/u.k,'TP')") % Planck temperature
    PlanckTemperature   = evalEverythere_("scd(u.TP,'PlanckTemperature')") 

    %---- pressure ----

    Pa                  = evalEverythere_("scd(u.N/u.sqm,'Pa')") % pascal
    pascal              = evalEverythere_("scd(u.Pa,'pascal')") 
    mPa                 = evalEverythere_("scd(1e-3*u.Pa,'mPa')") % millipascal
    millipascal         = evalEverythere_("scd(u.mPa,'millipascal')") 
    uPa                 = evalEverythere_("scd(1e-6*u.Pa,'uPa')") % micropascal
    micropascal         = evalEverythere_("scd(u.uPa,'micropascal')") 
    kPa                 = evalEverythere_("scd(1e3*u.Pa,'kPa')") % kilopascal
    kilopascal          = evalEverythere_("scd(u.kPa,'kilopascal')") 
    MPa                 = evalEverythere_("scd(1e6*u.Pa,'MPa')") % megapascal
    megapascal          = evalEverythere_("scd(u.MPa,'megapascal')") 
    GPa                 = evalEverythere_("scd(1e9*u.Pa,'GPa')") % gigapascal
    gigapascal          = evalEverythere_("scd(u.GPa,'gigapascal')") 
    bar                 = evalEverythere_("scd(1e5*u.Pa,'bar')") % bar
    mbar                = evalEverythere_("scd(1e-3*u.bar,'mbar')") % millibar
    millibar            = evalEverythere_("scd(u.mbar,'millibar')") 
    kbar                = evalEverythere_("scd(1e3*u.bar,'kbar')") % kilobar
    kilobar             = evalEverythere_("scd(u.kbar,'kilobar')") 
    atm                 = evalEverythere_("scd(101325*u.Pa,'atm')") % standard atmosphere
    atmosphere          = evalEverythere_("scd(u.atm,'atmosphere')") 
    standardAtmosphere  = evalEverythere_("scd(u.atm,'standardAtmosphere')") 
    at                  = evalEverythere_("scd(u.kgf/u.sqcm,'at')") % technical atmosphere
    technicalAtmosphere = evalEverythere_("scd(u.at,'technicalAtmosphere')") 
    torr                = evalEverythere_("scd(u.atm/760,'torr')") % torr
    Torr                = evalEverythere_("scd(u.torr,'Torr')") % torr
    mtorr               = evalEverythere_("scd(1e-3*u.torr,'mtorr')") % millitorr
    millitorr           = evalEverythere_("scd(u.mtorr,'millitorr')") 
    psi                 = evalEverythere_("scd(u.lbf/u.sqin,'psi')") % pound force per square inch
    poundPerSquareInch  = evalEverythere_("scd(u.psi,'poundPerSquareInch')") 
    ksi                 = evalEverythere_("scd(1e3*u.psi,'ksi')") % kip per square inch
    kipPerSquareInch    = evalEverythere_("scd(u.ksi,'kipPerSquareInch')") 
    Msi                 = evalEverythere_("scd(1e6*u.psi,'Msi')") % million pound force per square inch
    megapoundPerSquareInch = evalEverythere_("scd(u.Msi,'megapoundPerSquareInch')") 
    psf                 = evalEverythere_("scd(u.lbf/u.sqft,'psf')") % pound force per square foot
    poundPerSquareFoot  = evalEverythere_("scd(u.psf,'poundPerSquareFoot')") 
    ksf                 = evalEverythere_("scd(u.kip/u.sqft,'ksf')") % kip per square foot
    kipPerSquareFoot    = evalEverythere_("scd(u.ksf,'kipPerSquareFoot')") 
    Ba                  = evalEverythere_("scd(0.1*u.Pa,'Ba')") % barye
    barye               = evalEverythere_("scd(u.Ba,'barye')") 
    pz                  = evalEverythere_("scd(u.kPa,'pz')") % pièze
    pieze               = evalEverythere_("scd(u.pz,'pieze')") 
    mmHg                = evalEverythere_("scd(133.322387415*u.Pa,'mmHg')") % millimeter of mercury
    millimeterMercury   = evalEverythere_("scd(u.mmHg,'millimeterMercury')") 
    cmHg                = evalEverythere_("scd(10*u.mmHg,'cmHg')") % centimeter of mercury
    centimeterMercury   = evalEverythere_("scd(u.cmHg,'centimeterMercury')") 
    mHg                 = evalEverythere_("scd(1e3*u.mmHg,'mHg')") % meter of mercury
    meterMercury        = evalEverythere_("scd(u.mHg,'meterMercury')") 
    inHg                = evalEverythere_("scd(2.54*u.cmHg,'inHg')") % inch of mercury
    inchMercury         = evalEverythere_("scd(u.inHg,'inchMercury')") 
    ftHg                = evalEverythere_("scd(12*u.inHg,'ftHg')") % foot of mercury
    footMercury         = evalEverythere_("scd(u.ftHg,'footMercury')") 
    mmH20               = evalEverythere_("scd(u.kgf/u.sqm,'mmH20')") % millimeter of water (density 1 g/cc)
    mmAq                = evalEverythere_("scd(u.mmH20,'mmAq')") % millimeter of water
    millimeterWater     = evalEverythere_("scd(u.mmH20,'millimeterWater')") 
    cmH20               = evalEverythere_("scd(10*u.mmH20,'cmH20')") % centimeter of water
    cmAq                = evalEverythere_("scd(u.cmH20,'cmAq')") % centimeter of water
    centimeterWater     = evalEverythere_("scd(u.cmH20,'centimeterWater')") 
    mH20                = evalEverythere_("scd(1e3*u.mmH20,'mH20')") % meter of water
    mAq                 = evalEverythere_("scd(u.mH20,'mAq')") % meter of water
    meterWater          = evalEverythere_("scd(u.mH20,'meterWater')") 
    inH20               = evalEverythere_("scd(2.54*u.cmH20,'inH20')") % inch of water
    inAq                = evalEverythere_("scd(u.inH20,'inAq')") % inch of water
    inchWater           = evalEverythere_("scd(u.inH20,'inchWater')") 
    wc                  = evalEverythere_("scd(u.inH20,'wc')") % inch water column
    inchWaterColumn     = evalEverythere_("scd(u.wc,'inchWaterColumn')") 
    ftH20               = evalEverythere_("scd(12*u.inH20,'ftH20')") % foot of water
    ftAq                = evalEverythere_("scd(u.ftH20,'ftAq')") % foot of water
    footWater           = evalEverythere_("scd(u.ftH20,'footWater')") 

    %---- viscosity ----

    St                  = evalEverythere_("scd(u.sqcm/u.s,'St')") % stokes
    stokes              = evalEverythere_("scd(u.St,'stokes')") 
    cSt                 = evalEverythere_("scd(u.St/100,'cSt')") % centistokes
    centistokes       	= evalEverythere_("scd(u.cSt,'centistokes')") 
    newt                = evalEverythere_("scd(u.sqin/u.s,'newt')") % newt
    P                   = evalEverythere_("scd(u.Pa*u.s / 10,'P')") % poise
    poise               = evalEverythere_("scd(u.P,'poise')") 
    cP                  = evalEverythere_("scd(u.mPa*u.s,'cP')") % centipoise
    centipoise          = evalEverythere_("scd(u.cP,'centipoise')") 
    reyn                = evalEverythere_("scd(u.lbf*u.s/u.sqin,'reyn')") % reyn

    %---- power ----

    W                   = evalEverythere_("scd(u.J/u.s,'W')") % watt
    watt                = evalEverythere_("scd(u.W,'watt')") 
    kW                  = evalEverythere_("scd(1e3*u.W,'kW')") % kilowatt
    kilowatt            = evalEverythere_("scd(u.kW,'kilowatt')") 
    MW                  = evalEverythere_("scd(1e6*u.W,'MW')") % megawatt
    megawatt            = evalEverythere_("scd(u.MW,'megawatt')") 
    GW                  = evalEverythere_("scd(1e9*u.W,'GW')") % gigawatt
    gigawatt            = evalEverythere_("scd(u.GW,'gigawatt')") 
    TW                  = evalEverythere_("scd(1e12*u.W,'TW')") % terawatt
    terawatt            = evalEverythere_("scd(u.TW,'terawatt')") 
    mW                  = evalEverythere_("scd(1e-3*u.W,'mW')") % milliwatt
    milliwatt           = evalEverythere_("scd(u.mW,'milliwatt')") 
    uW                  = evalEverythere_("scd(1e-6*u.W,'uW')") % microwatt
    microwatt           = evalEverythere_("scd(u.uW,'microwatt')") 
    nW                  = evalEverythere_("scd(1e-9*u.W,'nW')") % nanowatt
    nanowatt            = evalEverythere_("scd(u.nW,'nanowatt')") 
    pW                  = evalEverythere_("scd(1e-12*u.W,'pW')") % picowatt
    picowatt            = evalEverythere_("scd(u.pW,'picowatt')") 
    hp                  = evalEverythere_("scd(550*u.ft*u.lbf/u.s,'hp')") % mechanical horsepower (550 ft-lbf/s)
    horsepower          = evalEverythere_("scd(u.hp,'horsepower')") 
    HP_I                = evalEverythere_("scd(u.hp,'HP_I')") % mechanical horsepower (550 ft-lbf/s)
    hpE                 = evalEverythere_("scd(746*u.W,'hpE')") % electrical horsepower
    HP_E                = evalEverythere_("scd(u.hpE,'HP_E')") % electrical horsepower
    electricalHorsepower = evalEverythere_("scd(u.hpE,'electricalHorsepower')") 
    PS                  = evalEverythere_("scd(75*u.kg*u.g0*u.m/u.s,'PS')") % metric horsepower (DIN 66036)
    HP                  = evalEverythere_("scd(u.PS,'HP')") % metric horsepower (DIN 66036)
    HP_DIN              = evalEverythere_("scd(u.PS,'HP_DIN')") % metric horsepower (DIN 66036)
    metricHorsepower    = evalEverythere_("scd(u.PS,'metricHorsepower')") 
    VA                  = evalEverythere_("scd(u.W,'VA')")
    volt_ampere         = evalEverythere_("scd(u.VA,'volt_ampere')")
    kVA                 = evalEverythere_("scd(1000*u.VA,'kVA')")
    kilovolt_ampere     = evalEverythere_("scd(u.kVA,'kilovolt_ampere')")
    var                 = evalEverythere_("scd(u.VA,'var')")
    volt_ampere_reactive = evalEverythere_("scd(u.var,'volt_ampere_reactive')")

    %---- current ----

    amp                 = evalEverythere_("scd(u.A,'amp')") % ampere
    ampere              = evalEverythere_("scd(u.A,'ampere')") 
    mA                  = evalEverythere_("scd(1e-3*u.A,'mA')") % milliampere
    milliampere         = evalEverythere_("scd(u.mA,'milliampere')") 
    uA                  = evalEverythere_("scd(1e-6*u.A,'uA')") % microampere
    microampere       	= evalEverythere_("scd(u.uA,'microampere')") 
    nA                  = evalEverythere_("scd(1e-9*u.A,'nA')") % nanoampere
    nanoampere          = evalEverythere_("scd(u.nA,'nanoampere')") 
    pA                  = evalEverythere_("scd(1e-12*u.A,'pA')") % picoampere
    picoampere          = evalEverythere_("scd(u.pA,'picoampere')") 
    kA                  = evalEverythere_("scd(1e3*u.A,'kA')") % kiloampere
    kiloampere          = evalEverythere_("scd(u.kA,'kiloampere')") 
    abA                 = evalEverythere_("scd(10*u.A,'abA')") % abampere
    abampere            = evalEverythere_("scd(u.abA,'abampere')") 
    Bi                  = evalEverythere_("scd(u.abA,'Bi')") % biot
    biot                = evalEverythere_("scd(u.Bi,'biot')") 

    %---- charge ----

    C                   = evalEverythere_("scd(u.A*u.s,'C')") % coulomb
    coulomb             = evalEverythere_("scd(u.C,'coulomb')") 
    mC                  = evalEverythere_("scd(1e-3*u.C,'mC')") % millicoulomb
    millicoulomb        = evalEverythere_("scd(u.mC,'millicoulomb')") 
    uC                  = evalEverythere_("scd(1e-6*u.C,'uC')") % microcoulomb
    microcoulomb        = evalEverythere_("scd(u.uC,'microcoulomb')") 
    nC                  = evalEverythere_("scd(1e-9*u.C,'nC')") % nanocoulomb
    nanocoulomb         = evalEverythere_("scd(u.nC,'nanocoulomb')") 
    pC                  = evalEverythere_("scd(1e-12*u.C,'pC')") % picocoulomb
    picocoulomb         = evalEverythere_("scd(u.pC,'picocoulomb')") 
    abC                 = evalEverythere_("scd(10*u.C,'abC')") % abcoulomb
    aC                  = evalEverythere_("scd(u.abC,'aC')") % abcoulomb
    abcoulomb           = evalEverythere_("scd(u.abC,'abcoulomb')") 
    statC               = evalEverythere_("scd(u.dyn^(1/2)*u.cm,'statC')") % statcoulomb
    statcoulomb         = evalEverythere_("scd(u.statC,'statcoulomb')") 
    Fr                  = evalEverythere_("scd(u.statC,'Fr')") % franklin
    franklin            = evalEverythere_("scd(u.Fr,'franklin')") 
    esu                 = evalEverythere_("scd(u.statC,'esu')") % electrostatic unit of charge
    electrostaticUnitCharge = evalEverythere_("scd(u.esu,'electrostaticUnitCharge')")
    mAh                 = evalEverythere_("scd(u.mA*u.hr,'mAh')") % milliamp-hour
    milliamp_hour       = evalEverythere_("scd(u.mAh,'milliamp_hour')")
    Ah                  = evalEverythere_("scd(u.A*u.hr,'Ah')") % amp-hour
    amp_hour            = evalEverythere_("scd(u.Ah,'amp_hour')") 

    %---- voltage ----

    V                   = evalEverythere_("scd(1*u.J/u.C,'V')") % volt
    volt                = evalEverythere_("scd(u.V,'volt')") 
    kV                  = evalEverythere_("scd(1e3*u.V,'kV')") % kilovolt
    kilovolt            = evalEverythere_("scd(u.kV,'kilovolt')") 
    MV                  = evalEverythere_("scd(1e6*u.V,'MV')") % megavolt
    megavolt            = evalEverythere_("scd(u.MV,'megavolt')") 
    GV                  = evalEverythere_("scd(1e9*u.V,'GV')") % gigavolt
    gigavolt            = evalEverythere_("scd(u.GV,'gigavolt')") 
    mV                  = evalEverythere_("scd(1e-3*u.V,'mV')") % millivolt
    millivolt           = evalEverythere_("scd(u.mV,'millivolt')") 
    uV                  = evalEverythere_("scd(1e-6*u.V,'uV')") % microvolt
    microvolt           = evalEverythere_("scd(u.uV,'microvolt')") 

    %---- resistance/conductance ----

    Ohm                 = evalEverythere_("scd(u.V/u.A,'Ohm')") % ohm
    GOhm                = evalEverythere_("scd(1e9*u.Ohm,'GOhm')") % gigaohm
    gigaohm            	= evalEverythere_("scd(u.GOhm,'gigaohm')") 
    MOhm                = evalEverythere_("scd(1e6*u.Ohm,'MOhm')") % megaohm
    megaohm             = evalEverythere_("scd(u.MOhm,'megaohm')") 
    kOhm                = evalEverythere_("scd(1e3*u.Ohm,'kOhm')") % kiloohm
    kiloohm           	= evalEverythere_("scd(u.kOhm,'kiloohm')") 
    mOhm                = evalEverythere_("scd(1e-3*u.Ohm,'mOhm')") % milliohm
    milliohm            = evalEverythere_("scd(u.mOhm,'milliohm')") 
    uOhm                = evalEverythere_("scd(1e-6*u.Ohm,'uOhm')") % microohm
    microohm            = evalEverythere_("scd(u.uOhm,'microohm')") 
    nOhm                = evalEverythere_("scd(1e-9*u.Ohm,'nOhm')") % nanoohm
    nanoohm             = evalEverythere_("scd(u.nOhm,'nanoohm')") 
    abOhm               = evalEverythere_("scd(u.nOhm,'abOhm')") % abohm
    Z0                  = evalEverythere_("scd(2*u.h_c*u.alpha/u.e^2,'Z0')") % characteristic impedance of vacuum
    impedanceOfVacuum  	= evalEverythere_("scd(u.Z0,'impedanceOfVacuum')") 
    R_K                 = evalEverythere_("scd(u.h_c/u.e^2,'R_K')") % von Klitzing constant
    vonKlitzingConstant = evalEverythere_("scd(u.R_K,'vonKlitzingConstant')") 
    R_K_90              = evalEverythere_("scd(25812.807*u.Ohm,'R_K_90')") % von Klitzing constant (conventional value)
    vonKlitzingConstant_conv = evalEverythere_("scd(u.R_K_90,'vonKlitzingConstant_conv')") 
    S                   = evalEverythere_("scd(1/u.Ohm,'S')") % siemens
    siemens             = evalEverythere_("scd(u.S,'siemens')") 
    mS                  = evalEverythere_("scd(1e-3*u.S,'mS')") % millisiemens
    millisiemens      	= evalEverythere_("scd(u.mS,'millisiemens')") 
    uS                  = evalEverythere_("scd(1e-6*u.S,'uS')") % microsiemens
    microsiemens        = evalEverythere_("scd(u.uS,'microsiemens')") 
    nS                  = evalEverythere_("scd(1e-9*u.S,'nS')") % nanosiemens
    nanosiemens         = evalEverythere_("scd(u.nS,'nanosiemens')") 
    G0                  = evalEverythere_("scd(2*u.e^2/u.h_c,'G0')") % conductance quantum
    conductanceQuantum  = evalEverythere_("scd(u.G0,'conductanceQuantum')") 

    %---- capacitance ----

    F                   = evalEverythere_("scd(u.A*u.s/u.V,'F')") % farad
    farad               = evalEverythere_("scd(u.F,'farad')") 
    mF                  = evalEverythere_("scd(1e-3*u.F,'mF')") % millifarad
    millifarad          = evalEverythere_("scd(u.mF,'millifarad')")
    uF                  = evalEverythere_("scd(1e-6*u.F,'uF')") % microfarad
    microfarad          = evalEverythere_("scd(u.uF,'microfarad')") 
    nF                  = evalEverythere_("scd(1e-9*u.F,'nF')") % nanofarad
    nanofarad           = evalEverythere_("scd(u.nF,'nanofarad')") 
    pF                  = evalEverythere_("scd(1e-12*u.F,'pF')") % picofarad
    picofarad           = evalEverythere_("scd(u.pF,'picofarad')")

    %---- inductance ----

    H                   = evalEverythere_("scd(u.Ohm*u.s,'H')") % henry
    henry               = evalEverythere_("scd(u.H,'henry')") 
    mH                  = evalEverythere_("scd(1e-3*u.H,'mH')") % millihenry
    millihenry          = evalEverythere_("scd(u.mH,'millihenry')") 
    uH                  = evalEverythere_("scd(1e-6*u.H,'uH')") % microhenry
    microhenry          = evalEverythere_("scd(u.uH,'microhenry')") 
    nH                  = evalEverythere_("scd(1e-9*u.H,'nH')") % nanohenry
    nanohenry           = evalEverythere_("scd(u.nH,'nanohenry')") 
    abH                 = evalEverythere_("scd(u.nH,'abH')") % abhenry
    abhenry             = evalEverythere_("scd(u.abH,'abhenry')") 
    kH                  = evalEverythere_("scd(1e3*u.H,'kH')") % kilohenry
    kilohenry           = evalEverythere_("scd(u.kH,'kilohenry')") 
    MH                  = evalEverythere_("scd(1e6*u.H,'MH')") % megahenry
    megahenry           = evalEverythere_("scd(u.MH,'megahenry')") 
    GH                  = evalEverythere_("scd(1e9*u.H,'GH')") % gigahenry
    gigahenry           = evalEverythere_("scd(u.GH,'gigahenry')") 

    %---- EM ----

    T                   = evalEverythere_("scd(1*u.N/(u.A*u.m),'T')") % tesla
    tesla               = evalEverythere_("scd(u.T,'tesla')") 
    Gs                  = evalEverythere_("scd(1e-4*u.T,'Gs')") % gauss
    gauss               = evalEverythere_("scd(u.Gs,'gauss')") 
    Wb                  = evalEverythere_("scd(u.V*u.s,'Wb')") % weber
    weber               = evalEverythere_("scd(u.Wb,'weber')") 
    Mx                  = evalEverythere_("scd(1e-8*u.Wb,'Mx')") % maxwell
    maxwell             = evalEverythere_("scd(u.Mx,'maxwell')") 
    mWb                 = evalEverythere_("scd(u.Wb/1000,'mWb')") % milliweber
    milliweber          = evalEverythere_("scd(u.mWb,'milliweber')") 
    uWb                 = evalEverythere_("scd(1e-6*u.Wb,'uWb')") % microweber
    microweber          = evalEverythere_("scd(u.uWb,'microweber')") 
    nWb                 = evalEverythere_("scd(1e-9*u.Wb,'nWb')") % nanoweber
    nanoweber           = evalEverythere_("scd(u.nWb,'nanoweber')") 
    Oe                  = evalEverythere_("scd(250/pi*u.A/u.m,'Oe')") % oersted
    oersted             = evalEverythere_("scd(u.Oe,'oersted')") 
    Gb                  = evalEverythere_("scd(2.5/pi*u.A,'Gb')") % gilbert
    gilbert             = evalEverythere_("scd(u.Gb,'gilbert')") 

    %---- non-dimensionals ----

    percent             = 0.01 % %
    pct                 = evalEverythere_("u.percent") % %
    permil              = 0.001 % ‰
    permill             = evalEverythere_("u.permil") % ‰
    permille            = evalEverythere_("u.permil") % ‰
    permyriad           = 1e-4 % permyriad
    bp                  = evalEverythere_("u.permyriad") % basis point
    basisPoint          = evalEverythere_("u.bp")
    ppm                 = 1e-6 % part per million
    partPerMillion      = evalEverythere_("u.ppm") 
    ppb                 = 1e-9 % part per billion
    partPerBillion      = evalEverythere_("u.ppb")
    ppt                 = 1e-12 % part per trillion
    partPerTrillion     = evalEverythere_("u.ppt")
    ppq                 = 1e-15 % part per quadrillion
    partPerQuadrillion  = evalEverythere_("u.ppq") 
    
    %---- angles ----
    % Note: angles are dimensionless

    rad                 = 1 % radian
    radian              = evalEverythere_("u.rad")
    sr                  = 1 % steradian
    steradian           = evalEverythere_("u.sr")
    turn                = 2*pi*evalEverythere_("u.rad") % turn
    rev                 = evalEverythere_("u.turn") % revolution = 2*pi radians
    revolution          = evalEverythere_("u.rev")
    deg                 = evalEverythere_("u.turn/360") % degree
    degree              = evalEverythere_("u.deg")
    arcmin              = evalEverythere_("u.deg/60") % arcminute
    arcminute           = evalEverythere_("u.arcmin")
    arcsec              = evalEverythere_("u.arcmin/60") % arcsecond
    arcsecond           = evalEverythere_("u.arcsec")
    grad                = evalEverythere_("u.turn/400") % gradian
    gradian             = evalEverythere_("u.grad")
    
    %---- rotational speed ----

    rpm                 = evalEverythere_("scd(u.rev/u.min,'rpm')") % revolution per minute
    revolutionPerMinute = evalEverythere_("scd(u.rpm,'revolutionPerMinute')") 
    rps                 = evalEverythere_("scd(u.rev/u.s,'rps')") % revolution per second
    revolutionPerSecond = evalEverythere_("scd(u.rps,'revolutionPerSecond')") 

    %---- velocity ----

    mps                 = evalEverythere_("scd(u.m/u.s,'mps')") % meter per second
    meterPerSecond      = evalEverythere_("scd(u.mps,'meterPerSecond')") 
    kyne                = evalEverythere_("scd(u.cm/u.s,'kyne')") % kyne
    Kyne                = evalEverythere_("scd(u.kyne,'Kyne')") % kyne
    fps                 = evalEverythere_("scd(u.ft/u.s,'fps')") % foot per second
    footPerSecond       = evalEverythere_("scd(u.fps,'footPerSecond')") 
    fpm                 = evalEverythere_("scd(u.ft/u.min,'fpm')") % foot per minute
    footPerMinute       = evalEverythere_("scd(u.fpm,'footPerMinute')") 
    kt                  = evalEverythere_("scd(u.nmi/u.hr,'kt')") % knot
    kn                  = evalEverythere_("scd(u.kt,'kn')") % knot
    kts                 = evalEverythere_("scd(u.kt,'kts')") % knot
    knot                = evalEverythere_("scd(u.kt,'knot')") 
    knot_UK             = evalEverythere_("scd(u.nm_UK/u.hr,'knot_UK')") % British imperial knot
    KTAS                = evalEverythere_("scd(u.kt,'KTAS')") % knot
    nmph                = evalEverythere_("scd(u.kt,'nmph')") % nautical mile per hour
    nauticalMilePerHour = evalEverythere_("scd(u.nmph,'nauticalMilePerHour')") 
    kph                 = evalEverythere_("scd(u.km/u.hr,'kph')") % kilometer per hour
    kmh                 = evalEverythere_("scd(u.kph,'kmh')") % kilometer per hour
    kps                 = evalEverythere_("scd(u.km/u.s,'kps')") % kilometer per second
    kilometerPerHour    = evalEverythere_("scd(u.kmh,'kilometerPerHour')") 
    mph                 = evalEverythere_("scd(u.mi/u.hr,'mph')") % mile per hour
    milePerHour         = evalEverythere_("scd(u.mph,'milePerHour')") 

    %---- volume flow rate ----

    cfm                 = evalEverythere_("scd(u.ft^3/u.min,'cfm')") % cubic foot per minute
    cubicFootPerMinute  = evalEverythere_("scd(u.cfm,'cubicFootPerMinute')") 
    cfs                 = evalEverythere_("scd(u.ft^3/u.s,'cfs')") % cubic foot per second
    cubicFootPerSecond  = evalEverythere_("scd(u.cfs,'cubicFootPerSecond')") 
    gpm                 = evalEverythere_("scd(u.gal/u.min,'gpm')") % US customary gallon per minute
    gallonPerMinute     = evalEverythere_("scd(u.gpm,'gallonPerMinute')") 
    gph                 = evalEverythere_("scd(u.gal/u.hr,'gph')") % US customary gallon per hour
    gallonPerHour       = evalEverythere_("scd(u.gph,'gallonPerHour')") 
    gpm_UK              = evalEverythere_("scd(u.gal_UK/u.min,'gpm_UK')") % British imperial gallon per minute
    lpm                 = evalEverythere_("scd(u.l/u.min,'lpm')") % liter per minute
    literPerMinute      = evalEverythere_("scd(u.lpm,'literPerMinute')") 

    %---- fuel economy ----

    l_100km             = evalEverythere_("scd(u.l/(100*u.km),'l_100km')") % liter per 100 km
    literPer100kilometer= evalEverythere_("scd(u.l_100km,'literPer100kilometer')") 
    mpg                 = evalEverythere_("scd(u.mi/u.gal,'mpg')") % mile per gallon
    milePerGallon       = evalEverythere_("scd(u.mpg,'milePerGallon')") 

    %---- Luminance etc. ----

    candela             = evalEverythere_("scd(u.cd,'candela')") 
    asb                 = evalEverythere_("scd(u.cd/u.sqm,'asb')") % apostilb
    apostilb            = evalEverythere_("scd(u.asb,'apostilb')") 
    sb                  = evalEverythere_("scd(u.cd/u.sqcm,'sb')") % stilb
    stilb               = evalEverythere_("scd(u.sb,'stilb')") 
    ph                  = evalEverythere_("scd(1e4*u.cd*u.sr/u.sqm,'ph')") % phot
    phot                = evalEverythere_("scd(u.ph,'phot')") 
    cp                  = evalEverythere_("scd(0.981*u.cd,'cp')") % candlepower
    candlepower         = evalEverythere_("scd(u.cp,'candlepower')") 
    lm                  = evalEverythere_("scd(u.cd*u.sr,'lm')") % lumen
    lumen               = evalEverythere_("scd(u.lm,'lumen')") 
    lx                  = evalEverythere_("scd(u.lm/u.sqm,'lx')") % lux
    lux                 = evalEverythere_("scd(u.lx,'lux')") 
    nx                  = evalEverythere_("scd(1e-3*u.lx,'nx')") % nox
    nox                 = evalEverythere_("scd(u.nx,'nox')") 

    %---- other derived SI ----

    mole                = evalEverythere_("scd(u.mol,'mole')") 
    kat                 = evalEverythere_("scd(u.mol/u.s,'kat')") % katal
    katal               = evalEverythere_("scd(u.kat,'katal')") 
    M                   = evalEverythere_("scd(u.mol/u.L,'M')") % molar
    molar               = evalEverythere_("scd(u.M,'molar')") 
    molarity           	= evalEverythere_("scd(u.M,'molarity')") % molarity
    Nms                 = evalEverythere_("scd(u.N*u.m*u.s,'Nms')") % newton-meter-second
    newton_meter_second = evalEverythere_("scd(u.Nms,'newton_meter_second')") 

    %---- radiation ----

    Gy                  = evalEverythere_("scd(u.J/u.kg,'Gy')") % gray
    gray                = evalEverythere_("scd(u.Gy,'gray')") 
    Sv                  = evalEverythere_("scd(u.J/u.kg,'Sv')") % sievert
    sievert             = evalEverythere_("scd(u.Sv,'sievert')") 
    Rad                 = evalEverythere_("scd(u.Gy/100,'Rad')") % absorbed radiation dose
    rem                 = evalEverythere_("scd(u.Sv/100,'rem')") % roentgen equivalent man
    roentgenEquivalentMan = evalEverythere_("scd(u.rem,'roentgenEquivalentMan')") 
    roentgen            = evalEverythere_("scd(2.58e-4*u.C/u.kg,'roentgen')") % roentgen
    Ly                  = evalEverythere_("scd(u.cal_th/u.sqcm,'Ly')") % langley
    lan                 = evalEverythere_("scd(u.Ly,'lan')") % langley
    langley            	= evalEverythere_("scd(u.lan,'langley')") 
    Bq                  = evalEverythere_("scd(1/u.s,'Bq')") % becquerel
    becquerel           = evalEverythere_("scd(u.Bq,'becquerel')") 
    Ci                  = evalEverythere_("scd(3.7e10*u.Bq,'Ci')") % curie
    curie               = evalEverythere_("scd(u.Ci,'curie')") 

    %---- constants ----
    
    i = 1i
    j = 1j
    pi = pi % Archimedes' constant ?
    tau = 2*pi
    phi = (1 + sqrt(5))/2 % golden ratio
    EulersNumber = exp(1) % ("e" is reserved for elementary charge)

    sigma_SB            = evalEverythere_("scd(pi^2/60*u.k^4/u.h_bar^3/u.c^2,'sigma_SB')") % Stefan-Boltzmann constant
    Stefan_BoltzmannConstant = evalEverythere_("scd(u.sigma_SB,'Stefan_BoltzmannConstant')") 
    mu_B                = evalEverythere_("scd(u.e*u.h_bar/(2*u.m_e),'mu_B')") % Bohr magneton
    BohrMagneton        = evalEverythere_("scd(u.mu_B,'BohrMagneton')") 
    mu_N                = evalEverythere_("scd(u.e*u.h_bar/(2*u.m_p),'mu_N')") % nuclear magneton
    nuclearMagneton     = evalEverythere_("scd(u.mu_N,'nuclearMagneton')") 
    ly                  = evalEverythere_("scd(u.c*u.year,'ly')") % light-year
    lightYear           = evalEverythere_("scd(u.ly,'lightYear')") % light-year
    mu0                 = evalEverythere_("scd(2*u.alpha*u.h_c/u.e^2/u.c,'mu0')") % vacuum permeability
    vacuumPermeability  = evalEverythere_("scd(u.mu0,'vacuumPermeability')") 
    eps0                = evalEverythere_("scd(u.c^-2/u.mu0,'eps0')") % vacuum permittivity
    vacuumPermittivity  = evalEverythere_("scd(u.eps0,'vacuumPermittivity')") 
    NAh                 = evalEverythere_("scd(u.N_A*u.h_c,'NAh')") % molar Planck constant
    molarPlanckConstant = evalEverythere_("scd(u.NAh,'molarPlanckConstant')") 
    M_u                 = evalEverythere_("scd(u.g/u.mol,'M_u')") % molar mass constant
    molarMassConstant   = evalEverythere_("scd(u.M_u,'molarMassConstant')") 
    K_J                 = evalEverythere_("scd(2*u.e/u.h_c,'K_J')") % Josephson constant
    JosephsonConstant   = evalEverythere_("scd(u.K_J,'JosephsonConstant')") 
    K_J_90              = evalEverythere_("scd(483597.9*u.Hz/u.V,'K_J_90')") % Josephson constant (conv. value)
    JosephsonConstant_conv = evalEverythere_("scd(u.K_J_90,'JosephsonConstant_conv')") 
    F_c                 = evalEverythere_("scd(u.e*u.N_A,'F_c')") % Faraday constant
    FaradayConstant     = evalEverythere_("scd(u.F_c,'FaradayConstant')") 
    c1                  = evalEverythere_("scd(2*pi*u.h_c*u.c^2,'c1')") % first radiation constant
    firstRadiationConstant = evalEverythere_("scd(u.c1,'firstRadiationConstant')")
    c2                  = evalEverythere_("scd(u.h_c*u.c/u.k,'c2')") % second radiation constant
    secondRadiationConstant = evalEverythere_("scd(u.c2,'secondRadiationConstant')") 
    b_prime             = evalEverythere_("scd(2.821439372122078893*u.k/u.h_c,'b_prime')") % Wien frequency displ. law const.
    WienFrequencyDisplacementLawConstant = evalEverythere_("scd(u.b_prime,'WienFrequencyDisplacementLawConstant')") 
    b_c                 = evalEverythere_("scd(u.h_c*u.c/(4.965114231744276303*u.k),'b_c')") % Wien wavelength displ. law const.
    WienWavelengthDisplacementLawConstant = evalEverythere_("scd(u.b_c,'WienWavelengthDisplacementLawConstant')") 
    R_air               = evalEverythere_("scd(287.05287*u.J/u.kg/u.K,'R_air')") % spec. gas const., air (ESDU 77022)
    specificGasConstant_air = evalEverythere_("scd(u.R_air,'specificGasConstant_air')") 
    R_bar               = evalEverythere_("scd(u.N_A*u.k,'R_bar')") % molar gas constant
    molarGasConstant    = evalEverythere_("scd(u.R_bar,'molarGasConstant')") 
    radarStatuteMile    = evalEverythere_("scd(2*u.mi/u.c,'radarStatuteMile')") 
    radarNauticalMile   = evalEverythere_("scd(2*u.NM/u.c,'radarNauticalMile')") 
    radarDataMile       = evalEverythere_("scd(2*u.dataMile/u.c,'radarDataMile')") 
    radarKilometer      = evalEverythere_("scd(2*u.km/u.c,'radarKilometer')") 

    %---- digital information ----

    nibble              = evalEverythere_("scd(4*u.bit,'nibble')") 
    B                   = evalEverythere_("scd(8*u.bit,'B')") % byte
    byte                = evalEverythere_("scd(u.B,'byte')") 
    octet               = evalEverythere_("scd(u.B,'octet')") % octet
    kB                  = evalEverythere_("scd(1e3*u.B,'kB')") % kilobyte
    kilobyte         	= evalEverythere_("scd(u.kB,'kilobyte')") 
    MB                  = evalEverythere_("scd(1e6*u.B,'MB')") % megabyte
    megabyte            = evalEverythere_("scd(u.MB,'megabyte')") 
    GB                  = evalEverythere_("scd(1e9*u.B,'GB')") % gigabyte
    gigabyte            = evalEverythere_("scd(u.GB,'gigabyte')") 
    TB                  = evalEverythere_("scd(1e12*u.B,'TB')") % terabyte
    terabyte            = evalEverythere_("scd(u.TB,'terabyte')") 
    PB                  = evalEverythere_("scd(1e15*u.B,'PB')") % petabyte
    petabyte            = evalEverythere_("scd(u.PB,'petabyte')") 
    EB                  = evalEverythere_("scd(1e18*u.B,'EB')") % exabyte
    exabyte             = evalEverythere_("scd(u.EB,'exabyte')") 
    Kibit               = evalEverythere_("scd(2^10*u.bit,'Kibit')") % kibibit
    kibibit             = evalEverythere_("scd(u.Kibit,'kibibit')") 
    KiB                 = evalEverythere_("scd(2^10*u.B,'KiB')") % kibibyte
    KB                  = evalEverythere_("scd(u.KiB,'KB')") % kibibyte
    kibibyte            = evalEverythere_("scd(u.KB,'kibibyte')") 
    Mibit               = evalEverythere_("scd(2^20*u.bit,'Mibit')") % mebibit
    mebibit             = evalEverythere_("scd(u.Mibit,'mebibit')") 
    MiB                 = evalEverythere_("scd(2^20*u.B,'MiB')") % mebibyte
    mebibyte            = evalEverythere_("scd(u.MiB,'mebibyte')") 
    Gibit               = evalEverythere_("scd(2^30*u.bit,'Gibit')") % gibibit
    gibibit             = evalEverythere_("scd(u.Gibit,'gibibit')") 
    GiB                 = evalEverythere_("scd(2^30*u.B,'GiB')") % gibibyte
    gibibyte            = evalEverythere_("scd(u.GiB,'gibibyte')") 
    Tibit               = evalEverythere_("scd(2^40*u.bit,'Tibit')") % tebibit
    tebibit             = evalEverythere_("scd(u.Tibit,'tebibit')") 
    TiB                 = evalEverythere_("scd(2^40*u.B,'TiB')") % tebibyte
    tebibyte            = evalEverythere_("scd(u.TiB,'tebibyte')") 
    Pibit               = evalEverythere_("scd(2^50*u.bit,'Pibit')") % pebibit
    pebibit             = evalEverythere_("scd(u.Pibit,'pebibit')") 
    PiB                 = evalEverythere_("scd(2^50*u.B,'PiB')") % pebibyte
    pebibyte            = evalEverythere_("scd(u.PiB,'pebibyte')") 
    Eibit               = evalEverythere_("scd(2^60*u.bit,'Eibit')") % exbibit
    exbibit             = evalEverythere_("scd(u.Eibit,'exbibit')") 
    EiB                 = evalEverythere_("scd(2^60*u.B,'EiB')") % exbibyte
    exbibyte            = evalEverythere_("scd(u.EiB,'exbibyte')") 
    bps                 = evalEverythere_("scd(u.bit/u.s,'bps')") % bit per second
    bitPerSecond        = evalEverythere_("scd(u.bps,'bitPerSecond')") 
    kbps                = evalEverythere_("scd(1e3*u.bps,'kbps')") % kilobit per second
    kilobitPerSecond   	= evalEverythere_("scd(u.kbps,'kilobitPerSecond')") 
    Mbps                = evalEverythere_("scd(1e6*u.bps,'Mbps')") % megabit per second
    megabitPerSecond    = evalEverythere_("scd(u.Mbps,'megabitPerSecond')") 
    Gbps                = evalEverythere_("scd(1e9*u.bps,'Gbps')") % gigabit per second
    gigabitPerSecond    = evalEverythere_("scd(u.Gbps,'gigabitPerSecond')") 
    Tbps                = evalEverythere_("scd(1e12*u.bps,'Tbps')") % terabit per second
    terabitPerSecond    = evalEverythere_("scd(u.Tbps,'terabitPerSecond')") 

    %---- currency ----
    % For display purposes - not for exchange rates.
    % See also mathworks.com/matlabcentral/fileexchange/47255

    cent                = evalEverythere_("scd(u.currency/100,'cent')") % cent (currency)
    Cent                = evalEverythere_("scd(u.cent,'Cent')") % cent (currency)
    pip                 = evalEverythere_("scd(u.cent/100,'pip')") % pip (currency)
    USD                 = evalEverythere_("scd(u.currency,'USD')") % currency
    EUR                 = evalEverythere_("scd(u.currency,'EUR')") % currency
    GBP                 = evalEverythere_("scd(u.currency,'GBP')") % currency
    JPY                 = evalEverythere_("scd(u.currency,'JPY')") % currency
    AUD                 = evalEverythere_("scd(u.currency,'AUD')") % currency
    CAD                 = evalEverythere_("scd(u.currency,'CAD')") % currency
    CHF                 = evalEverythere_("scd(u.currency,'CHF')") % currency
    CNY                 = evalEverythere_("scd(u.currency,'CNY')") % currency
    dollar              = evalEverythere_("scd(u.currency,'dollar')") % currency
    franc               = evalEverythere_("scd(u.currency,'franc')") % currency

    %---- used by Matlab's symunit but not here ----
    % gg - gauge
    % land - league
    % ha_US - US survey hectare
    % molecule
    % HP_UK - British imperial horsepower
    % PS_SAE - net horsepower (SAE J1349)
    % PS_DIN - horsepower (DIN 70020)
    % dry volumes
end

%% METHODS
methods
    %% Plotting and display:
    function disp(o)
        f = fieldnames(o);
        for iField = 1:length(f)
            thisField = u.(f{iField});
            if isa(thisField,'DimVar')
                thisField = scd(thisField);
            elseif isa(thisField,'OffsetDimVar')
                thisField = "value + " + string(thisField.offset);
            end
            uDisplayStruct.(f{iField}) = thisField;
        end
                
        try    
            dispdisp(uDisplayStruct);
            % mathworks.com/matlabcentral/fileexchange/48637
        catch
            builtin('disp',uDisplayStruct);
            
            url = 'http://www.mathworks.com/matlabcentral/fileexchange/48637/';
            dlCmd = sprintf('matlab:unzip(websave(tempname,''%s%s''),pwd);u',...
                url,'?download=true');
            
            warning('The function <a href="%s">%s</a> %s\n%s',...
                'www.mathworks.com/matlabcentral/fileexchange/48637',...
                'dispdisp',...
                'is recommended for display of physical units.',...
                ['<a href="' dlCmd ...
                '">Direct download of dispdisp into current directory</a>']);
        end
    end
end
end

%% Processing base units.
function U = buildCoreUnits(baseUnitSystem)
% import functions in case if repository has been includen in a package.
% if not - `import .*` does nothing 
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

coreBaseNames = {'m' 'kg' 's' 'A' 'K' 'mol' 'cd' 'bit' 'currency' 'unit'};

if ischar(baseUnitSystem) && strcmpi('none',baseUnitSystem)
    % Use normal variables - not DimVars - if baseUnitSystem is 'none'.
    U = cell2struct(num2cell(ones(size(coreBaseNames))),coreBaseNames,2);
    return
end

validateattributes(baseUnitSystem,{'cell'},{'size',[10,2]},...
    'u','baseUnitSystem');

if ~iscellstr(u.baseUnitSystem(:,1)) %#ok<ISCLSTR>
    error('First column of baseUnitSystem must be type char.')
end

baseValues = baseUnitSystem(:,2);
if ~all(cellfun('isclass', baseValues, 'double'))
    error('Second column of baseUnitSystem must contain doubles.')
end

expos = eye(numel(coreBaseNames));
for i = 1:numel(baseValues)
    U.(coreBaseNames{i}) = DimVar(expos(i,:),baseValues{i});
end
end

%% Local helper functions
function fun = getFullFunctionName_(fName)
%GETFULLFUNCTIONNAME_(fName) - resolves function name
%EXAMPLE:
%   %In case of this repository placed into package subfolder `./+somePkg`
%   fun = getFullFunctionName_("u.baseUnitSystem");
%   % fun = somePkg.u.baseUnitSystem

    sFullHandle = strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.');    %  e.g.: pkgName.subpkgName.fName
    if isempty(sFullHandle)
        fun = fName;
    else
        fun = sprintf('%s.%s',sFullHandle,fName);
    end
end

function res = evalEverythere_(func)
%EVALEVERYTHERE_(func) - resolves function name and calls it
%EXAMPLE:
%   %Evaluate full function name of `u.baseUnitSystem`
%   %call it wit arguments(`:,1`)
%   %transpose result
%   baseNames = evalEverythere_("u.baseUnitSystem(:,1)")';
arguments
    func    (1,:)   char
end
    fName = extractBefore(func,'(');
    if isempty(fName)
        % no arguments call like func
        fName = func;
        fArgs = '';
    else
        % call with arguments like func(arguments)
        [~,fArgs] = regexp(func,'(?<=\().*(?=\))','tokens','match');
        tUPath = getFullFunctionName_('u');
        tUPathShort = tUPath(1:end-2);
        if contains(fArgs,'OffsetDimVar') && ~isempty(tUPathShort)
            fArgs = strrep(fArgs,'OffsetDimVar',[tUPathShort '.OffsetDimVar']);
        end
        if isempty(fArgs)
            fArgs = '';
        elseif ~isempty(tUPath)
            fArgs = strsplit(fArgs{:},',');
            %parse arguments of the function recursively
            for i=1:numel(fArgs)
                if contains(fArgs{i},'u.')
                    fArgs{i} = strrep(fArgs{i},'u.',[tUPath '.']);
                end
            end
            %compose arguments string
            fArgs = strjoin(fArgs,',');
        end
    end
    fun = getFullFunctionName_(fName);
    if isempty(fArgs)
        fArgsString = '';
    else
        fArgsString = sprintf('(%s)',fArgs);
    end

    res = evalin('base',strrep(sprintf('%s%s',fun,fArgsString),':',''':'''));
end


%%
%   Original inspiration for this tool by Rob deCarvalho.
%     http://www.mathworks.com/matlabcentral/fileexchange/authors/22148
%     http://www.mathworks.com/matlabcentral/fileexchange/10070
