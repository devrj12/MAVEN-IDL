PRO SEA_reg3aa

; 3a edited for TJ -- moments calculated with mvn_sta_get_d1 resolution (16 sec)
;for plotting the data saved by SEA_reg1aa 


TIC

dir1 = '/home/joshi/Linux/Desktop/Research/IDL Commands/
;dir2 = '/home/joshi/Linux/Desktop/Research/IDL Commands/FILES_NR/
dir2 = '/home/joshi/Linux/Desktop/Research/IDL Commands/FILES2/

cd, dir2
RESTORE, '2016-08-20-3691a2.sav'
;RESTORE, '2017-06-20-5278a2.sav'
;RESTORE, '2015-08-26-1768a2.sav'
cd, dir1

startdate = time_double('2016-08-20/00:00')
ndays     =  1
timespan,startdate,ndays

A = time_double('2016-08-20/00:00:00')
B = time_double('2016-08-20/23:59:59')
Time = [A:B:1] 
	
mvn_swia_load_l2_data,/tplot,/loadall   ; , qlevel = 0.1
mvn_mag_load,'l2_1sec'
mvn_mag_load,'L2_FULL'
	
mk = mvn_spice_kernels(/load)
spice_vector_rotate_tplot,'mvn_B_1sec','MAVEN_MSO',check = 'MAVEN_SC_BUS' ; Creates tplot variable 'mvn_B_1sec_MAVEN_MSO'

;add position
spice_position_to_tplot,'MAVEN','MARS',frame = 'MSO'

orbdata = mvn_orbit_num()
store_data,'orbnum',orbdata.peri_time,orbdata.num,dlimit={ytitle:'Orbit'}
tplot,var_label='orbnum'

mvn_sta_l2_load
mvn_sta_l2_tplot

!p.background = 255
!p.color = 0

get_data,'mvn_swim_velocity_mso',data = VELOCITY_MSO ; for solar-wind velocity
tvectot,'mvn_swim_velocity_mso',newname='|vel_mso|'


Time_vmso = velocity_mso.x ; time-resolutions mostly 4 sec (but there are data gaps as well)

ind_v1 = (where(Time lt Time_vmso[0]))[-1] ; the index in Time where Time_vmso begins
ind_v2a = (where(Time_vmso gt Time[-1]) - 1); the index in Time_vmso where Time ends 
ind_v2 = ind_v2a[0]
ind_v2_ind = -[n_elements(ind_v2a) + 1]; the index in Time_vmso upto which Time elements are 'contained'

Time_vmso1 = round(Time_vmso[0:ind_v2]) ; it should be same as Time_vmso[0:ind_v2_ind]

ind0 = []
for el = 0,n_elements(Time_vmso1)-1 do begin 
      ind0 = [ind0, where(Time_vmso1[el] eq Time)]; indices in Time for equivalent elements in Time_vmso1 
endfor

get_data, 'mvn_B_1sec_MAVEN_MSO', data = Mag 
Mag_field = Mag.Y[ind0,*]  
 
Elec = crossproduct1(-VELOCITY_MSO.Y[0:ind_v2_ind,*], Mag_field) ; Time-resolution 4 sec (The unit of velocity_mso is km/sec and Mag_field
store_data, 'Elec',data = {x:VELOCITY_MSO.X[0:ind_v2_ind], y:Elec} ; is nanoTesla. So, the electric field here is in units of micro-volts/m.
tvectot,'Elec',newname='|Elec|' 
options,'|Elec|','yrange',[-2000,2000]
get_data,'Elec', data = Elec 

tvectot,'mvn_B_1sec_MAVEN_MSO',newname='|B|'
options,'|B|', 'yrange',[-10,20]

deno = sqrt(Mag.Y[*,0]^2 + Mag.Y[*,1]^2 + Mag.Y[*,2]^2)
numo = Mag.Y[*,0] 
;The angle (in degrees) whose cosine is the ratio(numo/deno) 
Angle = 180/!PI*ACOS(numo/deno)  
store_data,'Cone', data = {X:Mag.X,Y:Angle}
options,'Cone','ystyle',1  
                  			
get_data, 'orbnum', data = orbnum
Orbit = mvn_orbit_num(t = Time) 
Orbit2 = FLOOR(Orbit)
ORB_NO = Orbit2[UNIQ(Orbit2)]
options,'orbnum', 'yrange',[ORB_NO[0] - 10,ORB_NO[-1] + 10]

yline = fltarr(86400) + 15
yline1 = fltarr(86400) + 165
get_data,'Cone', data = Cone1, alim=lim                    
store_data,'line15',data={x:Cone1.x,y:yline},dlim=lim
store_data,'line135',data={x:Cone1.x,y:yline1},dlim=lim
options,'line15','colors','y' 
options,'line135','colors','y' 
store_data,'Cone_with_line', data = 'Cone line15 line135', dlim = lim
options,'Cone_with_line',ytitle,'Cone Angle'	

get_data, 'mvn_sta_d0_E', data = STA_CH0
get_data, 'mvn_sta_d1_E', data = STA_CH1

;to check adjacent orbit -- stop here !
STOP

	
FOREACH element, ORB_NO DO BEGIN

	IF element eq DEN_VEL.ORBIT THEN BEGIN	

		orb_ind = where(Orbit ge element and Orbit lt (element+1))
		Orbit1 = Orbit(orb_ind)
		Time1 = Time(orb_ind)
		
  		ind_o1 =  (where(Elec.X gt Time1[0]))[0] ; To recall Elec.X 'spans' Time_vmso1
		ind_o2 =  (where(Elec.X gt Time1[-1]))[0]-1 ; so these indices give where in Elec.X is Time1 contained
		store_data, 'Elec_o',data = {x:Elec.X[ind_o1:ind_o2], y:Elec.Y[ind_o1:ind_o2,*]}
		get_data,'Elec_o', data = Elec_o 
		
		ind01 = []
		for el1 = 0,n_elements(Elec_o.x)-1 do begin 
      			ind001 = where(Time1 eq round(Elec_o.x[el1]))
			if ind001 ne -1 then begin
				ind01 = [ind01, ind001];indices in Time1 (in this orbit) where Elec_o.x elements are 'contained'
			endif		
		endfor
		
		store_data,'Mag1', data = {x:Mag.x(orb_ind),y:Mag.y(orb_ind,*)}
		get_data,'Mag1',data = Mag1
		
		Time11 = Time1[ind01] ; Time1 in 4 sec resolution (solar-wind resolution)
		Mag11 = Mag1.Y[ind01,*]
		ind01a = where(round(Time_vmso) eq Time11[0]) 
		ind01b = where(round(Time_vmso) eq Time11[-1])
		store_data,'vel_swo', data = {x:velocity_mso.x(ind01a:ind01b),y:velocity_mso.y(ind01a:ind01b,*)}
		get_data,'vel_swo',data = vel_swo

		store_data,'Cone11',data = {x:Cone1.x(orb_ind),y:Cone1.y(orb_ind)}
		get_data,'Cone11',data = Cone11

		get_data,'mvn_swica_en_counts',data = Arch1
		get_data,'mvn_swics_en_counts',data = Sur1
		
		store_data,'Cone_b',data={x:Time1[Den_Vel.indall],y:Cone11.Y[Den_Vel.indall]},dlim=lim
		options,'Cone_b','colors','r' 
		options,'Cone_b','PSYM', 2 ;Asterisk
		options,'Cone_b','SYMSIZE',0.3

		;about elements satisfying radial conditions with the solar-wind region
		IF (n_elements(DEN_VEL.indall_cone) lt 2) THEN BEGIN	
			store_data,'Cone_with_line', data = 'Cone line15 line135 Cone_b', dlim = lim
		ENDIF ELSE BEGIN ; if radial conditions are satisfied			
			store_data,'Cone_c',data ={x:Time1[Den_Vel.indall_cone],y:Cone11.Y[Den_Vel.indall_cone]},dlim=dlim 
			options,'Cone_c','colors','g'
			options,'Cone_c','PSYM', 5 ;triangle
			options,'Cone_c','SYMSIZE',1.0
			store_data,'Cone_with_line', data = 'Cone line15 line135 Cone_b Cone_c', dlim = lim
		ENDELSE	
		
		IF ISA(Arch1,'STRUCT') THEN BEGIN

			Arch1.YTITLE = 'Energy A (ev)'
			store_data,'mvn_swica_en_counts',data = Arch1

			get_data,'mvn_swica_ph_counts',data = Arch2
			Arch2.YTITLE = 'Phi A'
			store_data,'mvn_swica_ph_counts', data = Arch2

			get_data,'mvn_swica_th_counts',data = Arch3
			Arch3.YTITLE = 'Theta A'
			store_data,'mvn_swica_th_counts', data = Arch3

			tplot ,['|B|','mvn_swica_en_counts','mvn_swica_ph_counts','mvn_swica_th_counts','mvn_swim_density','mvn_swim_velocity','mvn_sta_c6_E','mvn_sta_c6_M','mvn_sta_c0_H_E','mvn_sta_c0_L_E', 'Cone_with_line']
		ENDIF ELSE IF ISA(Sur1,'STRUCT') THEN BEGIN

			Sur1.YTITLE = 'Energy S (ev)'
			store_data,'mvn_swics_en_counts', data = Sur1

			get_data,'mvn_swics_ph_counts',data = Sur2
			Sur2.YTITLE = 'Phi S'
			store_data,'mvn_swics_ph_counts', data = Sur2

			get_data,'mvn_swics_th_counts',data = Sur3
			Sur3.YTITLE= 'Theta S'
			store_data,'mvn_swics_th_counts', data = Sur3	
		
			tplot ,['|B|','mvn_swics_en_counts','mvn_swics_ph_counts','mvn_swics_th_counts','mvn_swim_density','mvn_swim_velocity','mvn_sta_c6_E','mvn_sta_c6_M','mvn_sta_c0_H_E','mvn_sta_c0_L_E','Cone_with_line']
		ENDIF

		tlimit,([Time1[0],Time1[-1]])
		;tlimit,([Time1[12802], Time1[14602]])
		;tlimit,([Time1[1402], Time1[5602]])
		;([Time1[13582], Time1[13882]])

		TJ11 = Den_Vel.TJ11
		;O+
		vel_mso = spice_vector_rotate(DEN_VEL.VELO2,Time1[TJ11],'MAVEN_STATIC','MAVEN_MSO') ; this Time needs to be changed depending upon radial or non-radial conditions being checked ! (den_vel.velo2 is in km/sec and so should be vel_mso)
		vel_x = vel_mso[0,*]
		vel_y = vel_mso[1,*]
		vel_z = vel_mso[2,*]

		Magy = Mag1.Y ; magnetic field for this orbit in 1 sec resolution (in MSO coordinates !) 
		Mag1_x = (Magy[*,0])[TJ11]
		Mag1_y = (Magy[*,1])[TJ11]
		Mag1_z = (Magy[*,2])[TJ11]

		;Angle between velocity [for this orbit] and x in MSO : Blue
		Angle1 = []
		vel_mag1 = []
		for j = 0, (n_elements(vel_x) -1) do begin
			vel_mag = sqrt(vel_x[j]^2 + vel_y[j]^2 + vel_z[j]^2)  
			vel_mag1 = [vel_mag1,vel_mag]
			deno1 = vel_mag
			numo1 = vel_x[j]
			Angle1 = [Angle1, 180/!PI*ACOS(numo1/deno1)]
		endfor	

		;Angle between velocity [for this orbit] and local B : green
		Angle2 = []
		for j = 0, (n_elements(vel_x) -1) do begin
			vel_mag = sqrt(vel_x[j]^2 + vel_y[j]^2 + vel_z[j]^2) 
			Mag1_mag = sqrt(Mag1_x[j]^2 + Mag1_y[j]^2 + Mag1_z[j]^2)
			deno2 = vel_mag*Mag1_mag
			numo2 = vel_x[j]*Mag1_x[j]  + vel_y[j]*Mag1_y[j] + vel_z[j]*Mag1_z[j]
			Angle2 = [Angle2, 180/!PI*ACOS(numo2/deno2)]
		endfor

		;Angle between solar-wind velocity [for this orbit] and local B -- both in MSO
		Angle2a = []
		velsw_mag1 = []
		Mag11_mag1 = []
		velsw_x = vel_swo.y[*,0]
		velsw_y = vel_swo.y[*,1]
		velsw_z = vel_swo.y[*,2]
		Mag11_x = Mag11[*,0]
		Mag11_y = Mag11[*,1]
		Mag11_z = Mag11[*,2]
		for j = 0, (n_elements(Time11) - 1) do begin
		        velsw_mag = sqrt(velsw_x[j]^2 + velsw_y[j]^2 + velsw_z[j]^2) 
			velsw_mag1 = [velsw_mag1, velsw_mag]
			Mag11_mag = sqrt(Mag11_x[j]^2 + Mag11_y[j]^2 + Mag11_z[j]^2)
			Mag11_mag1 = [Mag11_mag1, Mag11_mag]
			deno2a = velsw_mag*Mag11_mag
			numo2a = velsw_x[j]*Mag11_x[j]  + velsw_y[j]*Mag11_y[j] + velsw_z[j]*Mag11_z[j]
			Angle2a = [Angle2a, 180/!PI*ACOS(numo2a/deno2a)]	
		endfor

		Angle2b = Angle2a*(!PI)/180
		g_rad = (32*1.66e-27*velsw_mag1*1e3*sin(Angle2b))/(1.6e-19*Mag11_mag1*1e-9)

		;ind01aa = (where(Time11 gt Time1[6230]))[0] ;change these indices 
 		;ind01ab = (where(Time11 lt Time1[7030]))[-1]
		
		;g_rad_orb = g_rad[ind01aa:ind01ab]

		;Angle between velocity [for this orbit] and Electric Field : red
		Elec_oy = Elec_o.Y ; Elec field for this orbit in the resolution of the Electric field data [4 sec] : make it 16 sec
		Elec_x = Elec_oy[*,0]; units -- micro-volts/m (km/sec multiplied by nanoTesla)
		Elec_y = Elec_oy[*,1]
		Elec_z = Elec_oy[*,2]

		Time1_tj1 = Time1[TJ11] ; times where moments have been calculated (in 16 sec resolution) in this orbit-time
		ind_t1 = (where(Time1_tj1 ge Time1[ind01[0]]))[0]
		ind_t2 = (where(Time1_tj1 le Time1[ind01[-1]]))[-1]
		Time1_tj2 = Time1_tj1[ind_t1:ind_t2] ; times (in 16 sec resolution) where Electric field data starts and ends in this orbit

		Elec_xt = interpol(Elec_x,Time1[ind01],Time1_tj2) ; 4 sec electric field data is interpolated to 16 sec resolution (16 
		Elec_yt = interpol(Elec_y,Time1[ind01],Time1_tj2) ; sec and 4 sec may not coincide)
		Elec_zt = interpol(Elec_z,Time1[ind01],Time1_tj2)
		
		; 1 sec resolution (Electric field starts and begins)
		ind_a = (where(Time1 ge Time1[ind01[0]]))[0]
		ind_b = (where(Time1 le Time1[ind01[-1]]))[-1]

		Elec_xt1 = interpol(Elec_x,Time1[ind01],Time1[ind_a:ind_b]) ; doesn't span the full orbit
		Elec_yt1 = interpol(Elec_y,Time1[ind01],Time1[ind_a:ind_b]) 
		Elec_zt1 = interpol(Elec_z,Time1[ind01],Time1[ind_a:ind_b])

		Elec_xt2 = TS_SMOOTH(Elec_xt1, 60)
		Elec_yt2 = TS_SMOOTH(Elec_yt1, 60)
		Elec_zt2 = TS_SMOOTH(Elec_zt1, 60)
		
		store_data, 'Elec_st', data = {x: Time1[ind_a:ind_b], y: [[Elec_xt2],[Elec_yt2],[Elec_zt2]]}
		tvectot, 'Elec_st', newname = '|Elec_st|'
		
		Mag1_xt = TS_SMOOTH(Mag1.Y[*,0],60) ; full orbit data
		Mag1_yt = TS_SMOOTH(Mag1.Y[*,1],60)
		Mag1_zt = TS_SMOOTH(Mag1.Y[*,2],60)

		store_data, 'Mag1_st', data = {x: Mag1.X, y: [[Mag1_xt],[Mag1_yt],[Mag1_zt]]}
		tvectot, 'Mag1_st', newname = '|Mag1_st|'

		vel0_x = vel_x[ind_t1:ind_t2] ; in km/sec
		vel0_y = vel_y[ind_t1:ind_t2] ; indices ind_t1 & ind_t2 : positions in Time1[TJ11] where the electric_field data begins
		vel0_z = vel_z[ind_t1:ind_t2] ; so these indices are 'covered' by Time1_tj2

		store_data, 'Den', data = {x:Time1[TJ11], y:den_vel.den2}
		store_data, 'vel', data = {x:Time1[TJ11], y:Transpose(vel_mso)} ;transpose to have the vector dimensions same as other vectors in MSO co-ordinates
		get_data,'vel', data = vel ; 16 sec resolution
		tvectot,'vel',newname='|Vel|'

		vel_xt = TS_SMOOTH(interpol(vel.Y[*,0],Time1[TJ11],Time1[ind_a:ind_b]),60)
		vel_yt = TS_SMOOTH(interpol(vel.Y[*,1],Time1[TJ11],Time1[ind_a:ind_b]),60)
		vel_zt = TS_SMOOTH(interpol(vel.Y[*,2],Time1[TJ11],Time1[ind_a:ind_b]),60) ; interpolated and smooth
		store_data, 'vel_st', data = {x: Time1[ind_a:ind_b], y: [[vel_xt],[vel_yt],[vel_zt]]}
		tvectot, 'vel_st', newname = '|vel_st|'
		get_data,'vel_st', data = vel_st ; 16 sec resolution

		store_data, 'Elec_velx', data = {x: Time1[ind_a:ind_b], y: [[Elec_xt2],[vel_xt]]}
		store_data, 'Elec_vely', data = {x: Time1[ind_a:ind_b], y: [[Elec_yt2],[vel_yt]]}
		store_data, 'Elec_velz', data = {x: Time1[ind_a:ind_b], y: [[Elec_zt2],[vel_zt]]}

		Angle3 = []
		Elec_mag1 = []
		for j = 0, (n_elements(vel0_x) -1) do begin
			vel0_mag = sqrt(vel0_x[j]^2 + vel0_y[j]^2 + vel0_z[j]^2) 
			Elec_mag = sqrt(Elec_xt[j]^2 + Elec_yt[j]^2 + Elec_zt[j]^2)
			Elec_mag1 = [Elec_mag1, Elec_mag]
			deno3 = vel0_mag*Elec_mag
			numo3 = vel0_x[j]*Elec_xt[j]  + vel0_y[j]*Elec_yt[j] + vel0_z[j]*Elec_zt[j]
			Angle3 = [Angle3, 180/!PI*ACOS(numo3/deno3)] ; this Angle3 shall be plotted against Elec_o.x
		endfor
		
		Angle1a = []
		for j = 0, (n_elements(vel0_x) -1) do begin
			vel0_maga = sqrt(vel0_x[j]^2 + vel0_y[j]^2 + vel0_z[j]^2) 
			deno1a = vel0_maga
			numo1a = vel0_x[j]
			Angle1a = [Angle1a, 180/!PI*ACOS(numo1a/deno1a)] ; this Angle1a is a test 
		endfor
		

		store_data,'Angle1',data  = {x:Time1[TJ11], y:Angle1} ; Time resolution 16 sec
		options,'Angle1','colors','b' 
		store_data,'Angle2',data  = {x:Time1[TJ11], y:Angle2} ; Time resolution 16 sec
		options,'Angle2','colors','g' 	
		store_data,'Angle3', data = {x:Time1_tj2,y:Angle3} ; Time resolution 16 sec
		options,'Angle3','colors','r' 
		store_data,'Angle1a',data  = {x:Time1_tj2, y:Angle1a}
		options, 'Angle1a','colors','k'
		store_data,'Angles', data = 'Angle1 Angle2 Angle3'
		options,'Angles', 'yrange',[0,180]
		options,'Angles','ytitle','O2+ velocity angle to :!c x in MSO (blue) !c B(green) !c Electric Field (red)'
		store_data, 'Den', data = {x:Time1[TJ11], y:den_vel.den2}
		store_data, 'vel', data = {x:Time1[TJ11], y:Transpose(vel_mso)} ; transpose to have the vector dimensions same as other vectors in MSO co-ordinates
		tvectot,'vel',newname='|Vel|'

		get_data,'MAVEN_POS_(MARS-MSO)',data = POS_MAVEN ; position
		p_d1 = abs(POS_MAVEN.X - Time1[0]) 
		p_d2 = where(p_d1 eq min(p_d1))
		p_d3 = abs(POS_MAVEN.X - Time1[-1])
		p_d4 = where(p_d3 eq min(p_d3))
		POS_MAVEN1_x = POS_MAVEN.X[p_d2[-1]:p_d4[0]]
		POS_MAVEN1_y = POS_MAVEN.Y[p_d2[-1]:p_d4[0],*]

		;get_data, 'mvn_B_1sec_MAVEN_MSO', data = Mag ; position
		m_d1 = abs(Mag.X - Time1[0]) 
		m_d2 = where(m_d1 eq min(m_d1))
		m_d3 = abs(Mag.X - Time1[-1])
		m_d4 = where(m_d3 eq min(m_d3))
		Mag_T_x = Mag.X[m_d2[-1]:m_d4[0]]
		Mag_T_y = Mag.Y[m_d2[-1]:m_d4[0],*]
		store_data, 'Mag_T', data = {x:Mag_T_x,y:Mag_T_y} ; magnetic data for this orbit

		sz_deno = sqrt(POS_MAVEN1_y[*,0]^2 + POS_MAVEN1_y[*,1]^2 + POS_MAVEN1_y[*,2]^2)
	        sz_numo = POS_MAVEN1_y[*,0]
		sz_Angle = 180/!PI*ACOS(sz_numo/sz_deno) 
		store_data, 'SZ_A', data = {x:POS_MAVEN1_x,y:sz_Angle}
		yline_sz = fltarr(n_elements(POS_MAVEN1_x)) + 90
		store_data,'line90',data={x:POS_MAVEN1_x,y:yline_sz},dlim=lim
		options,'line90','colors','r' 
		store_data,'SZ_A_with_line', data = 'SZ_A line90', dlim = lim

		store_data, 'ALT_MAVEN', data = {x:POS_MAVEN1_x,y:sz_deno - 3396.2}
		yline_alt = fltarr(n_elements(POS_MAVEN1_x)) + 350
		store_data,'line_a',data={x:POS_MAVEN1_x,y:yline_alt},dlim=lim
		options,'line_a','colors','r' 
		store_data,'ALT_MAVEN_line', data = 'ALT_MAVEN line_a', dlim = lim
		;options,'ALT_MAVEN_line', 'yrange',[0,500]
                
		;get_data,'|Vel|', data = vel_mag
		store_data,'Elec_z', data = {x: Time1[ind_a:ind_b], y: [Elec_zt2]}
		store_data,'Mag_yt', data = {x: Time1[ind_a:ind_b], y: [Mag1_yt[ind_a:ind_b]]}
		options,'Mag_yt','colors','b' 
		options,'Elec_z','colors','g' 
		
		store_data,'Elec_y', data = {x: Time1[ind_a:ind_b], y: [Elec_yt2]}
		store_data,'Mag_zt', data = {x: Time1[ind_a:ind_b], y: [Mag1_zt[ind_a:ind_b]]}
		options,'Mag_zt','colors','g' 
		options,'Elec_y','colors','b'  

		yline2 = fltarr(86400) + 0 
		store_data,'line0', data = {x:Time1[ind_a:ind_b],y:yline2},dlim=lim
 		options,'line0','colors','r' 
		store_data, 'Elec_linez', data = 'Elec_z line0', dlim = lim
		store_data, 'Mag_liney', data = 'Mag_yt line0', dlim = lim
		store_data, 'Elec_liney', data = 'Elec_y line0', dlim = lim
		store_data, 'Mag_linez', data = 'Mag_zt line0', dlim = lim
		store_data, 'MAVEN_line', data = 'MAVEN_POS_(MARS-MSO) line0', dlim = lim 

		;tplot,[132,4,5,6,13,14,44,45,[36,37],136,145]
		;tplot,[172,4,5,6,13,14,44,45,176,185]
		;tplot,[171,4,5,6,13,14,30,44,45,36,37,175]
		;tplot,[32,4,5,6,13,14,50,51,42,43,36,185]

		;tplot,[ '|B|',  'mvn_swica_en_counts', 'mvn_swica_ph_counts', 'mvn_swica_th_counts', 'mvn_swim_density',  'mvn_swim_velocity', 'mvn_sta_c6_E', 'mvn_sta_c6_M', 'mvn_sta_c0_H_E', 'mvn_sta_c0_L_E', 'Cone_with_line', 'Angles']
		;tplot,[ '|B|',  'mvn_swica_en_counts', 'mvn_swica_ph_counts', 'mvn_swica_th_counts', 'mvn_sta_c6_E', 'mvn_sta_c6_M', 'Cone_with_line', 'Angles', 'Den']

                 ; GENERALLY PLOTTED
		;tplot,[ '|B|',  'mvn_swica_en_counts', 'mvn_swica_ph_counts', 'mvn_swica_th_counts','mvn_sta_c6_E', 'mvn_sta_c6_M', 'mvn_sta_c0_H_E', 'mvn_sta_c0_L_E', 'Cone_with_line','Angles','Den','|Vel|','SZ_A_with_line', 'ALT_MAVEN_line','MAVEN_POS_(MARS-MSO)']
		;tplot,[ '|B|',  'mvn_swica_en_counts','mvn_sta_c6_M', 'Cone_with_line','Angles','Den','|Vel|','SZ_A_with_line', 'ALT_MAVEN_line','MAVEN_POS_(MARS-MSO)']
                 ;tplot,[ '|B|',  'mvn_swica_en_counts','mvn_sta_c6_M', 'Cone_with_line','Angles','Den','|Vel|','SZ_A_with_line', 'ALT_MAVEN_line','MAVEN_POS_(MARS-MSO)','mvn_sta_d1_D','mvn_sta_d1_A']

		;tplot,[ '|B|',  'mvn_swica_en_counts', 'mvn_swica_ph_counts', 'mvn_swica_th_counts','mvn_sta_c6_E', 'mvn_sta_c6_M', 'Angles','Cone_with_line']
	
		;tplot,[ '|B|',  'mvn_swica_en_counts', 'mvn_swica_ph_counts','mvn_swica_th_counts','mvn_swim_density','mvn_swim_velocity','mvn_sta_c6_M', 'Cone_with_line']
		;tplot,[ '|B|',  'mvn_swica_en_counts', 'mvn_sta_c6_M', 'Cone_with_line','Angles']
		;tplot,[ '|B|',  'mvn_swica_en_counts', 'mvn_sta_c6_M', 'Cone_with_line']
		

       ;smoothed;tplot,['|Mag1_st|', '|vel_st|','|Elec_st|','mvn_sta_c6_M','Elec_xt2a','vel_xt2a','Elec_yt2a','vel_yt2a', 'Elec_zt2a','vel_zt2a','MAVEN_POS_(MARS-MSO)'] 
		;tplot,['mvn_sta_c6_M','MAVEN_line','|Vel|','Elec_linez','Mag_liney', 'Elec_liney', 'Mag_linez','|Mag1_st|','|Elec_st|']
                ; same electric field in above two tplot commands as seen in 
				; store_data, 'Elec_zt2a', data = {x: Time1[ind_a:ind_b], y: Elec_zt2}
				; store_data, 'Elec_linez', data = 'Elec_z line0', dlim = lim
				; store_data,'Elec_z', data = {x: Time1[ind_a:ind_b], y: [Elec_zt2]}
		
	;non-smoothed;tplot,['mvn_sta_c6_M','mvn_swica_en_counts','mvn_sta_c6_E', 'mvn_sta_c0_H_E', 'mvn_sta_c0_L_E', 'MAVEN_line','Elec_linenz','Mag_lineny', 'Elec_lineny', 'Mag_linenz','|Mag1_st|','|Elec_st|','|Vel|'] 
		;tplot,['mvn_sta_c6_M', 'MAVEN_line','Elec_linenz','Mag_lineny', 'Elec_lineny', 'Mag_linenz','|Mag1_st|','|Elec_st|','|Vel|'] 
		;tplot,['|B|','|Vel|','|Elec|','mvn_sta_c6_M','Elec_ox','vel_x','Elec_oy','vel_y', 'Elec_oz','vel_z'] 

		;tplot,['mvn_sta_c6_M', 'MAVEN_line','Elec_linenz','Mag_lineny', 'Elec_lineny', 'Mag_linenz','|Mag1_st|','|Elec_st|','|Vel|','mvn_sta_d1_D', 'mvn_sta_d1_A']  
	

		STOP
	ENDIF	
	;STOP
ENDFOREACH

TOC

END

;TS = time_string(Time1[0],TFORMAT="YYYY-MM-DD-hh-mm-ss")
;TS1 = STRING(TS)
;TS3 = TS1 + 'a.png'
	
;path = '/home/joshi/Linux/Desktop/Research/Images/NEW/'		
;filename1 = path + TS3
;WRITE_PNG,filename1, TVRD(/TRUE)



