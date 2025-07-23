dt = monthday.(dts)
cdt=sort(unique(dt))
# climatological average call on monthday(dt)
function climprof(x, dt, cdt=sort(unique(dt)), bincond = (x,y) -> ==(x,y))
  ax = zeros(Float64, size(x,1), length(cdt))
  sx = zeros(Int64,   size(x,1), length(cdt))
  for I in CartesianIndices(x)
    i = I[1]; j = I[2]
    jc = bincond(dt[j], cdt)
    if any(jc) && !isnan(x[i,j]) && isfinite(x[i,j])
      ax[i,jc] += x[i,j]
      sx[i,jc] += 1
    end
  end
  return ax./sx
end
