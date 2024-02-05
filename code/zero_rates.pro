pro bond_template__define

  dimension = 12L
  template = {bond_template, $
    version:          0.0, $
    dataStart:        0L , $
    delimiter:        0B , $
    missingValue:     !values.f_nan, $
    commentSymbol:    '' , $
    fieldCount:       0L , $
    fieldTypes:       lonarr(dimension), $
    fieldNames:       strarr(dimension), $
    fieldLocations:   lonarr(dimension), $
    fieldGroups:      lonarr(dimension)  $
    }

end

pro zero_curve_template__define

  dimension = 2L
  template = {zero_curve_template, $
    version:          0.0, $
    dataStart:        0L , $
    delimiter:        0B , $
    missingValue:     !values.f_nan, $
    commentSymbol:    '' , $
    fieldCount:       0L , $
    fieldTypes:       lonarr(dimension), $
    fieldNames:       strarr(dimension), $
    fieldLocations:   lonarr(dimension), $
    fieldGroups:      lonarr(dimension)  $
    }

end

pro zero_surface_template__define

  dimension = 3L
  template = {zero_surface_template, $
    version:          0.0, $
    dataStart:        0L , $
    delimiter:        0B , $
    missingValue:     !values.f_nan, $
    commentSymbol:    '' , $
    fieldCount:       0L , $
    fieldTypes:       lonarr(dimension), $
    fieldNames:       strarr(dimension), $
    fieldLocations:   lonarr(dimension), $
    fieldGroups:      lonarr(dimension)  $
    }

end

pro zero_rates, time_min=time_min, time_max=time_max, T_min=T_min, T_max=T_max $
              , zero_min=zero_min, zero_max=zero_max, choice=choice, output=output

; Load the colors

loadcolors
black = 0
green = 4
blue  = 6
red   = 5

fieldcount_bond = 12L
fieldtypes_bond = [5L,5L,5L,5L,5L,5L,5L,5L,5L,5L,5L,5L]
fieldnames_bond = ['time', 'px_last','field3','field4','field5','bond_rate_last','field7','field8', $
                   'field9','field10','field11', 'zero_curve']
fieldlocations_bond = [0L,11L,19L,28L,37L,46L,55L,67L,79L,91L,103L,115L]
fieldgroups_bond = lindgen(fieldcount_bond)

bond_template = {bond_template, 1.0, 1L, 9B, !values.f_nan, '', fieldcount_bond, fieldtypes_bond, $
                fieldnames_bond, fieldlocations_bond, fieldgroups_bond}

bond_data = read_ascii('Bond_8.txt', template = bond_template)

fieldcount_zero_curve = 2L
fieldtypes_zero_curve = [5L,5L]
fieldnames_zero_curve = ['time_to_maturity', 'zero_rate']
fieldlocations_zero_curve = [0L,7L]
fieldgroups_zero_curve = lindgen(fieldcount_zero_curve)

zero_curve_template = {zero_curve_template, 1.0, 0L, 9B, !values.f_nan, '', fieldcount_zero_curve, $
fieldtypes_zero_curve, fieldnames_zero_curve, fieldlocations_zero_curve, fieldgroups_zero_curve}

zero_curve_data = read_ascii('zero_curve_8.txt', template = zero_curve_template)

fieldcount_zero_surface = 3L
fieldtypes_zero_surface = [5L,5L]
fieldnames_zero_surface = ['time','time_to_maturity', 'zero_rate']
fieldlocations_zero_surface = [0L,2L,9L]
fieldgroups_zero_surface = lindgen(fieldcount_zero_surface)

zero_surface_template = {zero_surface_template, 1.0, 0L, 9B, !values.f_nan, '', fieldcount_zero_surface, $
fieldtypes_zero_surface, fieldnames_zero_surface, fieldlocations_zero_surface, fieldgroups_zero_surface}

zero_surface_data = read_ascii('zero_surface_8.txt', template = zero_surface_template)

; Set the Postscript output

pson, filename=output, quiet=1B, page_size=[40.0,32.0]
device, /encapsulated, /cmyk

if (choice eq 'bond') then begin
    
    plot, bond_data.time, bond_data.bond_rate_last*100, thick=15.0 $
      , xtitle='Time (sec)', xstyle=1, xthick=15.0, xrange=[time_min,time_max] $
      , ytitle='Rates (%)', ystyle=1, ythick=15.0, yrange=[zero_min,zero_max] $
      , charthick=10.0, charsize=1.9, color=black, position=[0.125,0.135,0.88,0.98]
    
    xyouts, 0.73, 0.91, 'zero_rate', alignment=0.0, charthick=10, charsize=1.9, /normal
    plots, [0.77,0.85], [0.92,0.92], color=black, thick=15.0, /normal

    oplot, bond_data.time, bond_data.zero_curve*100, thick=15.0, color=blue
    xyouts, 0.74, 0.87, '!4q!3!N', alignment=0.0, charthick=10, charsize=1.9, /normal
    plots, [0.77,0.85], [0.88,0.88], color=blue, thick=15.0, /normal

endif

if (choice eq 'zero_curve') then begin
    
    plot, zero_curve_data.time, zero_curve.bond_rate_last*100, thick=15.0 $
      , xtitle='Time (sec)', xstyle=1, xthick=15.0, xrange=[time_min,time_max] $
      , ytitle='Rates (%)', ystyle=1, ythick=15.0, yrange=[zero_min,zero_max] $
      , charthick=10.0, charsize=1.9, color=black, position=[0.125,0.135,0.88,0.98]
    
    xyouts, 0.73, 0.91, 'zero_rate', alignment=0.0, charthick=10, charsize=1.9, /normal
    plots, [0.77,0.85], [0.92,0.92], color=black, thick=15.0, /normal

    oplot, bond_data.time, bond_data.zero_curve*100, thick=15.0, color=blue
    xyouts, 0.74, 0.87, '!4q!3!N', alignment=0.0, charthick=10, charsize=1.9, /normal
    plots, [0.77,0.85], [0.88,0.88], color=blue, thick=15.0, /normal

endif

; Close the Postscript device

psoff, quiet=0B

end

