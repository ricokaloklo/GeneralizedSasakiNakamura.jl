### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 1f0d7d7a-55a1-4a32-9a9f-2a4e0a30b101
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 5c4e0dc0-34b9-4e52-9c0f-d5f1f640a3b2
begin
    using GeneralizedSasakiNakamura
    using Printf
end

# ╔═╡ 4d30ec04-4ab2-4a89-9ff9-b6e2fd5119d0
HTML("""
<style>
pluto-output {
    max-width: 100% !important;
}
main {
    overflow-x: hidden !important;
}
</style>
""")

# ╔═╡ 728bf380-6bb2-4012-8e7e-97ce7c4d7b42
begin
    default_s = -2
    default_l = 2
    default_m = 2
    default_a = 0.7
    default_omega = 0.5
end

# ╔═╡ 0fbcfc49-81de-4452-87a1-5e358edc91f1
md""

# ╔═╡ 734a0d4b-d6bf-49c4-9498-e7255fd2e481
md""

# ╔═╡ f94a5856-1af9-414c-82ad-eedac5042fa6
md""

# ╔═╡ 19c21045-e740-4c9a-8625-b0bf419224ac
begin
    parse_int_control(x, default) = try
        parse(Int, string(x))
    catch
        default
    end

    parse_float_control(x, default) = try
        parse(Float64, string(x))
    catch
        default
    end

    function parse_demo_controls(raw)
        parts = split(string(raw), "|")
        length(parts) == 7 || return (default_s, default_l, default_m, default_a, default_omega, IN, 1600)
        s = parse_int_control(parts[1], default_s)
        l = parse_int_control(parts[2], default_l)
        m = clamp(parse_int_control(parts[3], default_m), -l, l)
        a = parse_float_control(parts[4], default_a)
        omega = parse_float_control(parts[5], default_omega)
        boundary = parts[6] == "UP" ? UP : IN
        ngrid = parse_int_control(parts[7], 1600)
        return (s, l, m, a, omega, boundary, ngrid)
    end

end

# ╔═╡ 45fb83ad-a8a9-4b35-8fc4-2f313d8f6bc7
md""

# ╔═╡ 3a61205f-8dcb-4073-aaf9-e3526d8e0973
md""

# ╔═╡ 894fc691-0060-4c6b-a06f-455dc374385c
md""

# ╔═╡ 1729d364-1888-4861-81cf-957fe3e783ec
md""

# ╔═╡ 2a59c0cf-618c-430c-818b-4b1e61cc69e0
begin
end

# ╔═╡ 69f946d2-9d5e-43e6-9d95-49ef54188436
md"""
"""

# ╔═╡ 93a65f2e-19b4-4a15-8210-70df148c519e
function radial_demo_data(; s, l, m, a, omega, boundary, ngrid)
    solve_time = @elapsed begin
        R = Teukolsky_radial(s, l, m, a, omega, boundary)
        X = GSN_radial(s, l, m, a, omega, boundary)
    end

    rsgrid = collect(range(-20, 100; length = ngrid))
    rgrid = r_from_rstar.(a, rsgrid)

    eval_time = @elapsed begin
        Rvals = R.(rgrid)
        Xvals = X.(rsgrid)
    end

    return (; R, X, rgrid, rsgrid, Rvals, Xvals, solve_time, eval_time)
end

# ╔═╡ f418c7fc-0e98-4a75-b01d-dced20f289f4
demo = nothing

# ╔═╡ 56c221a3-e0e7-47e4-bded-00fcfb8bc2d4
function svg_two_curves(xs, y1, y2; width = 900, height = 330, title = "", xlabel = "", ylabel = "")
    xmin, xmax = extrema(xs)
    ymin = min(minimum(y1), minimum(y2))
    ymax = max(maximum(y1), maximum(y2))
    ymin == ymax && (ymin -= 1; ymax += 1)
    pad_left = 70
    pad_right = 28
    pad_top = 58
    pad_bottom = 64
    plot_width = width - pad_left - pad_right
    plot_height = height - pad_top - pad_bottom

    x_to_px(x) = pad_left + plot_width * (x - xmin) / (xmax - xmin)
    y_to_px(y) = pad_top + plot_height * (ymax - y) / (ymax - ymin)
    points1 = join((
        @sprintf("%.3f,%.3f", x_to_px(xs[i]), y_to_px(y1[i]))
        for i in eachindex(xs)
    ), " ")
    points2 = join((
        @sprintf("%.3f,%.3f", x_to_px(xs[i]), y_to_px(y2[i]))
        for i in eachindex(xs)
    ), " ")

    zero_line = ymin < 0 < ymax ? @sprintf(
        "<line x1='%d' y1='%.3f' x2='%d' y2='%.3f' stroke='#cbd5e1' stroke-width='1' stroke-dasharray='5 4'/>",
        pad_left,
        y_to_px(0),
        width - pad_right,
        y_to_px(0),
    ) : ""
    xticks = join((
        begin
            x = xmin + (xmax - xmin) * k / 4
            px = x_to_px(x)
            @sprintf("<line x1='%.3f' y1='%d' x2='%.3f' y2='%d' stroke='#94a3b8' stroke-width='1'/><text x='%.3f' y='%d' text-anchor='middle' font-family='Helvetica, Arial, sans-serif' font-size='11' fill='#475569'>%.1f</text>", px, height - pad_bottom, px, height - pad_bottom + 6, px, height - pad_bottom + 22, x)
        end
        for k in 0:4
    ), "")
    yticks = join((
        begin
            y = ymin + (ymax - ymin) * k / 4
            py = y_to_px(y)
            @sprintf("<line x1='%d' y1='%.3f' x2='%d' y2='%.3f' stroke='#e2e8f0' stroke-width='1'/><text x='%d' y='%.3f' text-anchor='end' dominant-baseline='middle' font-family='Helvetica, Arial, sans-serif' font-size='11' fill='#475569'>%.2g</text>", pad_left, py, width - pad_right, py, pad_left - 8, py, y)
        end
        for k in 0:4
    ), "")

    return """
    <svg width="$width" height="$height" viewBox="0 0 $width $height" style="display:block;width:100%;height:auto;" xmlns="http://www.w3.org/2000/svg">
      <rect x="0" y="0" width="$width" height="$height" fill="#fbfbf8"/>
      <text x="$(width / 2)" y="28" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="20" font-weight="700" fill="#172033">$title</text>
      $yticks
      $xticks
      $zero_line
      <line x1="$pad_left" y1="$(height - pad_bottom)" x2="$(width - pad_right)" y2="$(height - pad_bottom)" stroke="#334155" stroke-width="1.2"/>
      <line x1="$pad_left" y1="$pad_top" x2="$pad_left" y2="$(height - pad_bottom)" stroke="#334155" stroke-width="1.2"/>
      <text x="$(width / 2)" y="$(height - 16)" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#334155">$xlabel</text>
      <text x="20" y="$(height / 2)" transform="rotate(-90 20 $(height / 2))" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#334155">$ylabel</text>
      <polyline points="$points1" fill="none" stroke="#2667ff" stroke-width="2.2"/>
      <polyline points="$points2" fill="none" stroke="#d97706" stroke-width="2.2"/>
      <rect x="$(width - 170)" y="46" width="130" height="54" rx="8" fill="rgba(255,255,255,0.82)" stroke="#cbd5e1"/>
      <line x1="$(width - 154)" y1="64" x2="$(width - 116)" y2="64" stroke="#2667ff" stroke-width="2.5"/>
      <text x="$(width - 106)" y="68" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#172033">real</text>
      <line x1="$(width - 154)" y1="84" x2="$(width - 116)" y2="84" stroke="#d97706" stroke-width="2.5"/>
      <text x="$(width - 106)" y="88" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#172033">imag</text>
    </svg>
    """
end

# ╔═╡ 755b4134-b6ed-4e81-8011-7f3d45e2dc54
function format_complex_box(label, z)
    @sprintf(
        "<div class='amp-box'><div class='amp-label'>%s amplitude</div><div class='amp-value'>%.6e %+.6e i</div></div>",
        label,
        real(z),
        imag(z),
    )
end

# ╔═╡ 4f683a59-1104-43ad-9175-4dd15a4525a4
function demo_card(demo; s, l, m, a, omega, boundary)
    r_svg = svg_two_curves(
        demo.rgrid,
        real.(demo.Rvals),
        imag.(demo.Rvals);
        title = "Teukolsky function",
        xlabel = "r / M",
        ylabel = "R(r)",
    )
    x_svg = svg_two_curves(
        demo.rsgrid,
        real.(demo.Xvals),
        imag.(demo.Xvals);
        title = "GSN function",
        xlabel = "r* / M",
        ylabel = "X(r*)",
    )
    mode_text = @sprintf("s = %d, l = %d, m = %d, a = %.2f, omega = %.2f, boundary = %s", s, l, m, a, omega, string(boundary))
    solve_time_text = @sprintf("Solve equation on-the-fly: %.3f ms", 1000 * demo.solve_time)
    eval_time_text = @sprintf("Evaluation on %d grid points: %.3f ms", length(demo.rgrid), 1000 * demo.eval_time)
    return HTML("""
    <style>
      .gsn-demo-card {
        width: min(930px, calc(100vw - 48px));
        max-width: 100%;
        box-sizing: border-box;
        border: 1px solid #d6d3ca;
        border-top: 0;
        border-radius: 0 0 18px 18px;
        padding: 14px 14px 16px 14px;
        background: linear-gradient(180deg, #fffefb 0%, #f5f1e8 100%);
        box-shadow: 0 16px 42px rgba(15, 23, 42, 0.12);
        font-family: Helvetica, Arial, sans-serif;
      }
      .gsn-mode-line {
        display: flex;
        justify-content: space-between;
        gap: 12px;
        margin: 4px 8px 12px 8px;
        color: #243447;
        font-size: 12px;
        line-height: 1.2;
      }
      .gsn-mode-text {
        white-space: nowrap;
      }
      .gsn-time {
        color: #52616f;
        font-size: 11px;
        line-height: 1.2;
        white-space: nowrap;
        text-align: right;
      }
      .plot-pair {
        display: grid;
        grid-template-columns: 1fr;
        gap: 10px;
      }
      .amp-grid {
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 10px;
        margin: 12px 8px 0 8px;
      }
      .amp-box {
        border: 1px solid #d6d3ca;
        border-radius: 12px;
        background: rgba(255,255,255,0.76);
        padding: 10px 12px;
      }
      .amp-label {
        font-size: 12px;
        color: #64748b;
        margin-bottom: 6px;
      }
      .amp-value {
        font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
        font-size: 12px;
        color: #172033;
        overflow-wrap: anywhere;
      }
    </style>
    <div class='gsn-demo-card'>
      <div class='gsn-mode-line'>
        <div class='gsn-mode-text'>$mode_text</div>
        <div class='gsn-time'><div>$solve_time_text</div><div>$eval_time_text</div></div>
      </div>
      <div class='plot-pair'>$r_svg $x_svg</div>
      <div class='amp-grid'>
        $(format_complex_box("incidence", demo.R.incidence_amplitude))
        $(format_complex_box("reflection", demo.R.reflection_amplitude))
        $(format_complex_box("transmission", demo.R.transmission_amplitude))
      </div>
    </div>
    """)
end

# ╔═╡ 7e7b5749-b83c-4439-8770-c9b64360e4c9
@bind controls_raw HTML("""
<div class='gsn-bond-root'>
  <style>
  .gsn-bond-root {
      width: min(930px, calc(100vw - 48px));
      max-width: 100%;
      box-sizing: border-box;
      margin: 8px 0 0 0;
      font-family: Helvetica, Arial, sans-serif;
  }
  .pluto-controls {
      box-sizing: border-box;
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px 16px;
      align-items: center;
      padding: 12px 16px;
      border-radius: 18px 18px 0 0;
      background: #202020;
      color: #d6d6d6;
      font-weight: 700;
  }
  .pluto-controls label {
      display: grid;
      grid-template-columns: auto 1fr;
      gap: 8px;
      align-items: center;
  }
  .pluto-controls input[type='range'] {
      width: 100%;
      accent-color: #f7f7f7;
  }
  .pluto-controls select {
      width: 100%;
  }
  @media (max-width: 720px) {
      .pluto-controls {
          grid-template-columns: repeat(2, minmax(0, 1fr));
      }
  }
  </style>
  <div class='pluto-controls'>
    <label>s <select data-name='s'><option value='-2' selected>-2</option><option value='-1'>-1</option><option value='0'>0</option><option value='1'>1</option><option value='2'>2</option></select></label>
    <label>l <input data-name='l' type='range' min='2' max='30' step='1' value='2'></label>
    <label>m <input data-name='m' type='range' min='-2' max='2' step='1' value='2'></label>
    <label>a <input data-name='a' type='range' min='0.0' max='0.99' step='0.01' value='0.7'></label>
    <label>omega <input data-name='omega' type='range' min='0.05' max='1.0' step='0.01' value='0.5'></label>
    <label>boundary <select data-name='boundary'><option value='IN' selected>IN</option><option value='UP'>UP</option></select></label>
    <label>grid <input data-name='grid' type='range' min='400' max='3000' step='100' value='1600'></label>
  </div>
  <script>
  (() => {
      const root = currentScript.closest(".gsn-bond-root")
      const cell = root.closest("pluto-cell")
      if (cell) {
          const input = cell.querySelector("pluto-input")
          const runarea = cell.querySelector("pluto-runarea")
          const shoulder = cell.querySelector("pluto-shoulder")
          if (input) input.style.display = "none"
          if (runarea) runarea.style.display = "none"
          if (shoulder) shoulder.style.display = "none"
      }
      const panel = root.querySelector(".pluto-controls")
      const read = name => panel.querySelector("[data-name='" + name + "']").value
      const lSlider = panel.querySelector("[data-name='l']")
      const mSlider = panel.querySelector("[data-name='m']")
      const syncMRange = () => {
          const l = parseInt(lSlider.value, 10)
          const m = Math.max(-l, Math.min(l, parseInt(mSlider.value, 10)))
          mSlider.min = String(-l)
          mSlider.max = String(l)
          mSlider.value = String(m)
      }
      const centerCard = () => {
          const card = document.querySelector(".gsn-demo-card")
          if (card) {
              card.scrollIntoView({ block: "center", inline: "nearest" })
          }
      }
      const update = keepCentered => {
          syncMRange()
          root.value = [read("s"), read("l"), read("m"), read("a"), read("omega"), read("boundary"), read("grid")].join("|")
          root.dispatchEvent(new CustomEvent("input", { bubbles: true }))
          root.dispatchEvent(new Event("change", { bubbles: true }))
          if (keepCentered) {
              setTimeout(centerCard, 80)
              setTimeout(centerCard, 350)
              setTimeout(centerCard, 900)
          }
      }
      panel.querySelectorAll("input, select").forEach(el => {
          el.addEventListener("input", () => update(true))
          el.addEventListener("change", () => update(true))
      })
      syncMRange()
      update(false)
  })()
  </script>
</div>
""")

# ╔═╡ 4888058e-51a4-42b4-8d6d-7de35af7a601
begin
    current_s, current_l, current_m, current_a, current_omega, current_boundary, current_ngrid = parse_demo_controls(controls_raw)
    current_demo = radial_demo_data(;
        s = current_s,
        l = current_l,
        m = current_m,
        a = current_a,
        omega = current_omega,
        boundary = current_boundary,
        ngrid = current_ngrid,
    )

    card = demo_card(
        current_demo;
        s = current_s,
        l = current_l,
        m = current_m,
        a = current_a,
        omega = current_omega,
        boundary = current_boundary,
    )
    card
end

# ╔═╡ 07cc551e-a67a-4dd5-9177-913b3a352a23
md""

# ╔═╡ Cell order:
# ╠═1f0d7d7a-55a1-4a32-9a9f-2a4e0a30b101
# ╠═5c4e0dc0-34b9-4e52-9c0f-d5f1f640a3b2
# ╟─4d30ec04-4ab2-4a89-9ff9-b6e2fd5119d0
# ╠═728bf380-6bb2-4012-8e7e-97ce7c4d7b42
# ╟─0fbcfc49-81de-4452-87a1-5e358edc91f1
# ╟─734a0d4b-d6bf-49c4-9498-e7255fd2e481
# ╟─f94a5856-1af9-414c-82ad-eedac5042fa6
# ╠═19c21045-e740-4c9a-8625-b0bf419224ac
# ╟─45fb83ad-a8a9-4b35-8fc4-2f313d8f6bc7
# ╟─3a61205f-8dcb-4073-aaf9-e3526d8e0973
# ╟─894fc691-0060-4c6b-a06f-455dc374385c
# ╟─1729d364-1888-4861-81cf-957fe3e783ec
# ╠═2a59c0cf-618c-430c-818b-4b1e61cc69e0
# ╟─69f946d2-9d5e-43e6-9d95-49ef54188436
# ╠═93a65f2e-19b4-4a15-8210-70df148c519e
# ╠═f418c7fc-0e98-4a75-b01d-dced20f289f4
# ╠═56c221a3-e0e7-47e4-bded-00fcfb8bc2d4
# ╠═755b4134-b6ed-4e81-8011-7f3d45e2dc54
# ╠═4f683a59-1104-43ad-9175-4dd15a4525a4
# ╠═7e7b5749-b83c-4439-8770-c9b64360e4c9
# ╠═4888058e-51a4-42b4-8d6d-7de35af7a601
# ╟─07cc551e-a67a-4dd5-9177-913b3a352a23
