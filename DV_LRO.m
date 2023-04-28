clear all
clc

%% LOI DV computation

date = [2009 6 23 11 27 00]; % Arrival date at moon
day = date2mjd2000(date);
[r_moon, v_moon] = ephMoon(day); %position and velocity vector of moon at arrival of s/c
r2 = r_moon;
vers_r2 = r2./norm(r2); % unit vector of moon position vector
vers_r1 = -vers_r2;% start of transfer orbit from a position with opposite unit vector wrt to moon at arrival date
r1 = vers_r1.*(185+astroConstants(23));
mu = astroConstants(13);
mu_m = astroConstants(20);

% HOHMANN TRANSFER
a_transfer = (norm(r2)+norm(r1))/2; %semi-major axis of transfer orbit
e_transfer = (norm(r2)-norm(r1))/(norm(r2)+norm(r1)); % eccentricity of transfer orbit
h_transfer = sqrt(mu*a_transfer*(1-e_transfer^2)); % angular momentum of transfer orbit
v_apo_transfer = (mu/h_transfer)*(1-e_transfer); %velocity of s/c at apocentre (moon position)


v_moon_vers = v_moon./norm(v_moon);
v_rel = v_apo_transfer - norm(v_moon); % relative velocity between s/c and moon
v_r = norm(v_rel);
RM = astroConstants(30);


%LOI 1  capture (216x3081.4 km altitude) at pericentre (216km altitude)
a = (RM+216 + RM+3081.4)/2;
e = ((RM+3081.4) - (RM+216))/ ((RM+3081.4) + (RM+216));
h = sqrt(mu_m*a*(1-e^2));
v_ins = (mu_m/h)*(1+e); %velocity at pericentre of 216x3081.4 km
v_peri_hyp = sqrt(v_r^2 + (2*mu_m)/(RM+216)); %velocity at pericentre of hyperbola
DV = abs(v_peri_hyp-v_ins);  


%LOI 3    (216x740 km altitude) with manoeuvre at pericentre (216km altitude)
a_3 = (RM+216 + RM+740)/2;
e_3 = ((RM+740) - (RM+216))/ ((RM+740) + (RM+216));
h_3 = sqrt(mu_m*a_3*(1-e_3^2));
v_peri_3 = (mu_m/h_3)*(1+e_3); % velocity at pericentre of 216x740 km
DV3 = abs(v_ins - v_peri_3);

%LOI 4 circularization at 216 km altitude
v_circ = sqrt(mu_m/(RM+216));
DV4 = abs(v_peri_3 - v_circ);

%LOI5 commissioning orbit (216x30 km altitude) with manoeuvre at pericentre
% of 216x740 km, that is the apocentre of the commissioning orbit
a_c = (RM+216 + RM+30)/2;
e_c = ((RM+216) - (RM+30))/ ((RM+216) + (RM+30));
h_c = sqrt(mu_m*a_c*(1-e_c^2));
v_apo_comm = (mu_m/h_c)*(1-e_c); % velocity at apocentre of commissioning orbit
DV5 = abs(v_circ- v_apo_comm);

DV_tot = DV+DV3+DV4+DV5; % TOTAL COMPUTED DV FOR THE LOI MANOEUVRES

%% MOI (mission orbit insertion, considered as 50km of altitude circular orbit)

v_peri_comm = (mu_m/h_c)*(1+e_c); %pericentre velocity of commissioning orbit
a_m_insertion = (RM+50 + RM+30)/2;
e_m_insertion = ((RM+50) - (RM+30))/ ((RM+30) + (RM+50));
h_m = sqrt(mu_m*a_m_insertion*(1-e_m_insertion^2)); 
v_peri_mission = (mu_m/h_m)*(1+e_m_insertion); % velocity at pericentre of mission orbit
v_apo_mission = (mu_m/h_m)*(1-e_m_insertion); % velocity at apocentre of mission orbit
DV_mission1 = abs(v_peri_mission-v_peri_comm); %To become 30x50 km altitude manoeuvre at pericentre of commissioning orbit

v_circ_m = sqrt(mu_m/(RM+50)); % velocity in circular orbit
DV_mission2 = abs(v_apo_mission- v_circ_m); 

DV_mission = DV_mission2+ DV_mission1; % TOTAL DV OF MOI MANOEUVRES







