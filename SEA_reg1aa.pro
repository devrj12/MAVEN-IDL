PRO SEA_reg1aa
TIC
;This is for all orbits (which have some portion in solar-wind region) [but not only which have radial portion] in it! That's the edit from SEA_reg1.pro. And it will be run by SEA_reg3aa.pro
;SEA_reg1a after subtraction (of both high-counts and low-counts and also including TJ[every 16th sec moment calculations]).

FOR nj = 0,1 DO BEGIN
	
	startdate = time_string(time_double('2015-08-17/00:00')  + nj*86400)
	ndays     =  1
	timespan,startdate,ndays

	A = time_double('2015-08-17/00:00:00') + nj*86400
	B = time_double('2015-08-17/23:59:59') + nj*86400
	
	mvn_swia_load_l2_data,/tplot,/loadall
	mvn_mag_load,'l2_1sec'
	mvn_mag_load,'L2_FULL'
	
	mk = mvn_spice_kernels(/load)
	spice_vector_rotate_tplot,'mvn_B_1sec','MAVEN_MSO',check = 'MAVEN_SC_BUS'

	;add position
	spice_position_to_tplot,'MAVEN','MARS',frame = 'MSO'

	orbdata = mvn_orbit_num()
	store_data,'orbnum',orbdata.peri_time,orbdata.num,dlimit={ytitle:'Orbit'}
	tplot,var_label='orbnum'

	;mvn_sta_l2_load
	;mvn_sta_l2_tplot
	mvn_sta_l2_load,iv_level = 3
	mvn_sta_l2_tplot,/replace
	
	!p.background = 255
	!p.color = 0

	tvectot,'mvn_B_1sec_MAVEN_MSO',newname='|B|'
	options,'|B|', 'yrange',[-10,20]

	get_data, 'mvn_B_1sec_MAVEN_MSO', data = Mag 
	deno = sqrt(Mag.Y[*,0]^2 + Mag.Y[*,1]^2 + Mag.Y[*,2]^2)
	numo = Mag.Y[*,0] 
	Angle = 180/!PI*ACOS(numo/deno)  
	store_data,'Cone', data = {X:Mag.X,Y:Angle}
	options,'Cone','ystyle',1  
	
	get_data, 'orbnum', data = orbnum

	Time = [A:B:1] 
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
	
	STOP
	
	FOREACH element, ORB_NO DO BEGIN
	
	IF element eq 1719 then begin; if only a chosen orbit

		orb_ind = where(Orbit ge element and Orbit lt (element+1))
		Orbit1 = Orbit(orb_ind)
		Time1 = Time(orb_ind)
		
		store_data,'Cone11',data = {x:Cone1.x(orb_ind),y:Cone1.y(orb_ind)}
		get_data,'Cone11',data = Cone11

		get_data,'mvn_swica_en_counts',data = Arch1
		get_data,'mvn_swics_en_counts',data = Sur1
		get_data,'mvn_swim_density',data = Swim1 
		IF ISA(Swim1,'STRUCT') THEN BEGIN
			Time_swim = Swim1.X 
			w = where(Time_swim ge Time1[0] and Time_swim le Time1[-1],nel)
		ENDIF
		
		IF ISA(Arch1,'STRUCT') THEN BEGIN

			Arch1.YTITLE = 'Energy A (ev)'
			store_data,'mvn_swica_en_counts',data = Arch1

			get_data,'mvn_swica_ph_counts',data = Arch2
			Arch2.YTITLE = 'Phi A'
			store_data,'mvn_swica_ph_counts', data = Arch2

			get_data,'mvn_swica_th_counts',data = Arch3
			Arch3.YTITLE = 'Theta A'
			store_data,'mvn_swica_th_counts', data = Arch3

			tplot ,['|B|','mvn_swica_en_counts','mvn_swica_ph_counts','mvn_swica_th_counts','mvn_swim_density','mvn_swim_velocity','orbnum','mvn_sta_c6_E','mvn_sta_c6_M','mvn_sta_c0_H_E','mvn_sta_c0_L_E', 'Cone_with_line']
		ENDIF ELSE IF ISA(Sur1,'STRUCT') THEN BEGIN

			Sur1.YTITLE = 'Energy S (ev)'
			store_data,'mvn_swics_en_counts', data = Sur1

			get_data,'mvn_swics_ph_counts',data = Sur2
			Sur2.YTITLE = 'Phi S'
			store_data,'mvn_swics_ph_counts', data = Sur2

			get_data,'mvn_swics_th_counts',data = Sur3
			Sur3.YTITLE= 'Theta S'
			store_data,'mvn_swics_th_counts', data = Sur3	
		
			tplot ,['|B|','mvn_swics_en_counts','mvn_swics_ph_counts','mvn_swics_th_counts','mvn_swim_density','mvn_swim_velocity','orbnum','mvn_sta_c6_E','mvn_sta_c6_M','mvn_sta_c0_H_E','mvn_sta_c0_L_E','Cone_with_line']
		ENDIF
		
		tlimit,([Time1[0],Time1[-1]])
		IF [(ISA(Arch1,'STRUCT') eq 1) or (ISA(Sur1,'STRUCT') eq 1)] and [(typename(nel) ne 'UNDEFINED')] THEN BEGIN ; if data exists
			IF (nel gt 1) THEN BEGIN
			mvn_swia_regid
			get_data, 'regid', data = regid 
			regidy1 = regid.Y[*,0]
			reg_in  = where(regidy1 eq 1) ;solar wind region
			ENDIF
		ENDIF
		STOP
		; if there are elements in the solar-wind region
		IF [ISA(reg_in) eq 1] and [n_elements(reg_in) GE 2] THEN BEGIN
			regidy11 = regidy1(reg_in)
			n1 = n_elements(regidy1)*1.0
			n2 = n_elements(regidy11)*1.0
			n_fraction1 = n2/n1
			
			;STOP
			regidx1 = regid.X ; this is time from the "new region"
			regidx11 = regidx1(reg_in) ; solar-wind region in the new region
			regidx_f = floor(regidx11)
			reg_diff = -ts_diff(regidx11,1) ; difference between elements
			reg_diff1 = reg_diff[0:-2]; omitting last element as it is zero	
		
			; getting 'correct' elements
			regidx_f1 = floor(regidx11[0])
			FOR del = 0,n_elements(reg_diff1)-1  DO BEGIN
				regidx_f1 = [regidx_f1, regidx_f1[-1]+round(reg_diff1[del])]
			ENDFOR
			; getting indices in Time1 for the (solar-wind part in the new region)
			JJ1 = 0
			WHILE (JJ1 GE 0) DO BEGIN
				;ind1 = where(Time1 eq floor(regidx11[JJ1])) 
				ind1 = where(Time1 eq regidx_f1[JJ1]) 
				IF (ind1 ne -1) THEN BREAK
				JJ1 = JJ1 + 1
			ENDWHILE
	
			JJ2 = -1
			WHILE (JJ2 LE -1) DO BEGIN
				;ind2 = where(Time1 eq floor(regidx11[JJ2]))
				ind2 = where(Time1 eq regidx_f1[JJ2]) 
				IF (ind2 ne -1) THEN BREAK
				JJ2 = JJ2 - 1
			ENDWHILE
                
			;get all indices in Time1 : start with the first index and take every 4th index         
			indall = indgen((ind2-ind1)/4 + 1,start=ind1,increment=4); indices in Time1 where the new region (with solar-wind) starts and ends
			Time1_indall = floor(Time1(indall))
			;compare Time1_indall and regidx_f : setsubtract gives an array of elements which are uncommon in Time1_indall, regidx_f	
			SS = setsubtract(Time1_indall,regidx_f1)	
			IF (n_elements(SS) eq 1) and (SS[0] eq -1) THEN BEGIN	
				Cone_11 = Cone11.Y[indall] ; indices wise - Cone11.Y is equivalent to that of Time1.
				cone_ind = where(Cone_11 le 15 or Cone_11 ge 165)
				indall_cone = indall[cone_ind] ; indices which satisfy radial conditions (in solar wind)
				Cone_13 = Cone_11[cone_ind] ; values which satisfy radial conditions (in solar wind) (or, Cone11.Y[indall_cone])
				Time1_cone = Time1[indall_cone]; times which satisfy radial conditions (in solar wind)
				store_data,'Cone_b',data={x:Time1_indall,y:Cone_11},dlim=lim
			ENDIF ELSE BEGIN
				new_array = []
				FOREACH element2, regidx_f1 DO BEGIN
					diff = Time1_indall - element2
					diff1 = abs(diff)
					ind_ss = where(diff1 eq min(diff1))
					new_array = [new_array, ind_ss] ; 
				ENDFOREACH
				indall = indall[new_array]
				Cone_11 =  Cone11.Y(indall) ; Cone values in solar-wind region.Cone_11 should have same n_elements as regidx_f1. That's the check.
				cone_ind = where(Cone_11 le 15 or Cone_11 ge 165)
				;indall_conen = indall[new_array]; indices for being in solar wind region
				indall_cone = indall[cone_ind];indices [for this orbit] which satisfy radial conditions (in solar wind)
				Cone_13 = Cone_11[cone_ind] ; values which satisfy radial conditions (in solar wind)
				Time1_cone = Time1[indall_cone]; times which satisfy radial conditions (in solar wind)

				store_data,'Cone_b',data={x:Time1_indall[new_array],y:Cone_11},dlim=lim	
			ENDELSE

			n_1 = n_elements(Cone_11)*1.0	; solar wind region elements
			n_2 = n_elements(cone_ind)*1.0  ; radial condition region elements

			store_data,'Cone_a',data={x:Time1,y:Cone11.Y},dlim=lim
			options,'Cone_a','colors','b'
			options,'Cone_a','PSYM', 1 ; Plus sign(+)
			options,'Cone_a','SYMSIZE',0.2 
			options,'Cone_b','colors','r' 
			options,'Cone_b','PSYM', 2 ;Asterisk
			options,'Cone_b','SYMSIZE',0.3
			
			;about elements satisfying radial conditions with the solar-wind region
			IF (n_elements(cone_ind) eq 1) and (cone_ind[0] eq -1) THEN BEGIN	
				n_2 = 0			
				n_fraction2 = n_2/n_1
				store_data,'Cone_with_line', data = 'Cone line15 line135 Cone_b', dlim = lim
			ENDIF ELSE BEGIN ; if radial conditions are satisfied
				n_fraction2 = n_2/n_1
				store_data,'Cone_c',data ={x:Time1_cone,y:Cone_13},dlim=dlim 
				options,'Cone_c','colors','g'
				options,'Cone_c','PSYM', 5 ;triangle
				options,'Cone_c','SYMSIZE',1.0
				store_data,'Cone_with_line', data = 'Cone line15 line135 Cone_b Cone_c', dlim = lim
			ENDELSE
			
			den1 = []
			velo1 = []
			den2 = []
			velo2 = []
			den3 = []
			velo3 = []

			get_data, 'mvn_sta_d0_E', data = STA_CH0 ; d0 has 128 sec resolution
			get_data, 'mvn_sta_d1_E', data = STA_CH1 ; d1 has 16 sec resolution
			
			STOP
			IF isa(sta_ch1,'struct') eq 0 THEN continue ; if sta_ch1 isn't a structure, then go to next iteration
			sta_time = STA_CH1.X

			;IF ISA(sta_ch1,'STRUCT') THEN BEGIN 
			;	sta_time = STA_CH1.X
			;ENDIF ELSE IF ISA(sta_ch0,'STRUCT') THEN BEGIN
			;	sta_time = STA_CH0.X
                        ;ENDIF 

			TJ11 = []
			FOR tj = 0,n_elements(Time1)-1 DO BEGIN 
				print,tj
				if tj mod 16 ne 0 then continue
				TSD = sta_time - Time1[tj]
				TSD_wa = where(abs(TSD) eq min(abs(TSD)))
 				TSD_wb = abs(sta_time[TSD_wa] - Time1[tj])

				if TSD_wb gt 16 then continue
				
				dat  = mvn_sta_get_d1(Time1[tj]) 			
				mass = dat.mass_arr

				TJ11 = [TJ11, tj]

				;ind_tj1 = where(dat.DATA[*,*,0] gt (6270/8))
				
				;if total(ind_tj1) ne -1 then begin
				;	print, 'yes'
				;	for tj1 = 1,7 do begin
				;		AD1 = dat.DATA[*,*,tj1]
				;		AD2 = AD1[ind_tj1]
				;		AD2 = 0
				;		AD1[ind_tj1] = AD2
				;		dat.DATA[*,*,tj1] = AD1
						
				;		ADc1 = dat.CNTS[*,*,tj1]
				;		ADc2 = ADc1[ind_tj1]
				;		ADc2 = 0
				;		ADc1[ind_tj1] = ADc2
				;		dat.CNTS[*,*,tj1] = ADc1
				;	endfor
				;endif

				;ind_tj2 = where(dat.DATA[*,*,0] gt 0 and dat.DATA[*,*,0] lt (6270/8))
				;if total(ind_tj2) ne -1 then begin
					;print, 'yes_a'
					;for tj2 = 1,7 do begin
                                        ;	DD1 = dat.DATA[*,*,tj2]
					;	DD2 = DD1[ind_tj2]
					;	DD2a = dat.DATA[*,*,0]
					;	DD2b = DD2a[ind_tj2]
					;	DD2 = DD2 - 0.02*DD2b
					;	ind_tj3 = where(DD2 lt 0)
					;	DD2[ind_tj3] = 0
					;	DD1[ind_tj2] = DD2
					;	dat.DATA[*,*,tj2] = DD1

                                        ;	DDc1 = dat.CNTS[*,*,tj2]
					;	DDc2 = DDc1[ind_tj2]
					;	DDc2a = dat.CNTS[*,*,0]
					;	DDc2b = DDc2a[ind_tj2]
					;	DDc2 = DDc2 - 0.02*DDc2b
					;	ind_tjc3 = where(DDc2 lt 0)
					;	DDc2[ind_tjc3] = 0
					;	DDc1[ind_tj2] = DDc2
					;	dat.CNTS[*,*,tj2] = DDc1
					;endfor
				;endif

				;if tj eq 2853 then begin
				;	STOP					
				;endif

				;O+
				den1  = [den1, n_4d(dat,MASS=[12, 20],ENERGY = [25, 20000])] ; row is okay for it
				velo1 = [[velo1], [v_4d(dat,MASS=[12,20],ENERGY = [25, 20000])]] ; such bracketing to make them in columns
				
				;O2+
				den2  = [den2, n_4d(dat,MASS=[28,36],ENERGY = [25, 20000])]
				velo2 = [[velo2], [v_4d(dat,MASS=[28,36],ENERGY = [25, 20000])]] 

				;CO2+
				den3  = [den3, n_4d(dat,MASS=[40,48],ENERGY = [25, 20000])]
				velo3 = [[velo3], [v_4d(dat,MASS=[40,48],ENERGY = [25, 20000])]] 
			ENDFOR
			STOP
			IF TJ11 ne !NULL THEN BEGIN
				Den_Vel = Create_Struct('Date', STRMID(startdate,0,10), 'Orbit', element, 'TJ11', TJ11, 'Time1', Time1, 'Time1_cone', Time1_cone, 'orb_ind', orb_ind,'indall', indall, 'indall_cone', indall_cone, 'Den1', den1, 'Velo1', velo1, 'Den2', den2, 'velo2', velo2, 'Den3', den3, 'velo3', velo3)
				;Den_Vel = Create_Struct('Date', STRMID(startdate,0,10), 'Orbit', element, 'Time1', Time1, 'Time1_cone', Time1_cone, 'orb_ind', orb_ind,'indall', indall, 'indall_cone', indall_cone, 'Den1', den1, 'Velo1', velo1, 'Den2', den2, 'velo2', velo2, 'Den3', den3, 'velo3', velo3)
				
				AB = STRMID(startdate,0,10) + '-' + STRTRIM(element,2) + 'a2.sav'
				AC = '/home/joshi/Linux/Desktop/Research/IDL Commands/FILES2/'+ AB  
				IF FILE_TEST(AC) eq 0 THEN BEGIN					
					SAVE, FILENAME = AC, Den_Vel
				ENDIF
			ENDIF ; TJ11 check

		tplot
		STOP
		ENDIF ; reg_in check (if there are times in solar-wind in this orbit)
		reg_in = !NULL
		Arch1 = !NULL
		Sur1 = !NULL
		Swim1 = !NULL
		nel = !NULL	
	ENDIF	; for chosen orbit
	ENDFOREACH ; for each orbit
	;ENDIF
	;STOP
ENDFOR ; for each day
TOC
END

;rd = spice_vector_rotate(DEN_VEL.VELO1,time_string(DEN_VEL.TIME),'MAVEN_STATIC','MAVEN_MSO')
