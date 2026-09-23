! 本模块汇集表面发射相关的计算、随机采样与输入解析工具，供初始化、时间推进和后处理调用。
! 模块默认隐藏内部实现，只公开下方列出的状态码及过程，以减少调用方对实现细节的依赖。
module surface_emission_utils
   use, intrinsic :: iso_fortran_env, only: real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mt19937_64, only: genrand64_real3
   implicit none
   private

   integer, parameter, public :: SEEM_OK = 0
   integer, parameter, public :: SEEM_OUT_OF_RANGE = 1
   integer, parameter, public :: SEEM_ERR_INVALID = 2
   integer, parameter, public :: SEEM_ERR_OVERFLOW = 3

   ! SEEM_OUT_OF_RANGE 表示插值能量超出数据表范围；此时仍返回最近端点的产额值。
   ! 其余非零状态表示输入无效或计算超出可表示范围，调用方应检查 ierr 后再使用输出值。
   public :: accumulate_surface_emission
   public :: interpolate_see_yield
   public :: sample_poisson_count
   public :: compute_see_macro_mean
   public :: compute_kinetic_energy_ev
   public :: compute_surface_current_density
   public :: parse_surface_current_definition
   public :: parse_secondary_emission_definition
   public :: read_see_yield_table

contains

   ! 将观察窗口内进入计算域的总电荷换算为面电流密度，保留输入电荷的正负号。
   ! area 和 window_seconds 必须为正；非有限输入返回无效状态，中间结果溢出则返回溢出状态。
   subroutine compute_surface_current_density(charge_into_domain, area, window_seconds, current_density, ierr)
      real(real64), intent(in) :: charge_into_domain, area, window_seconds
      real(real64), intent(out) :: current_density
      integer, intent(out) :: ierr

      current_density = 0.0_real64
      ierr = SEEM_ERR_INVALID
      if (.not. ieee_is_finite(charge_into_domain) .or. .not. ieee_is_finite(area) .or. &
          .not. ieee_is_finite(window_seconds)) return
      if (area <= 0.0_real64 .or. window_seconds <= 0.0_real64) return
      current_density = charge_into_domain / area
      if (.not. ieee_is_finite(current_density)) then
         current_density = 0.0_real64
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      current_density = current_density / window_seconds
      if (.not. ieee_is_finite(current_density)) then
         current_density = 0.0_real64
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      ierr = SEEM_OK
   end subroutine compute_surface_current_density

   ! 用经典动能 1/2*m*v^2 计算单粒子动能，并除以元电荷换算成电子伏特。
   ! 质量和元电荷必须为正，速度分量及其他实数输入必须是有限值。
   subroutine compute_kinetic_energy_ev(mass, velocity, elementary_charge, energy_ev, ierr)
      real(real64), intent(in) :: mass, velocity(3), elementary_charge
      real(real64), intent(out) :: energy_ev
      integer, intent(out) :: ierr

      energy_ev = 0.0_real64
      ierr = SEEM_ERR_INVALID
      if (.not. ieee_is_finite(mass) .or. .not. all(ieee_is_finite(velocity)) .or. &
          .not. ieee_is_finite(elementary_charge)) return
      if (mass <= 0.0_real64 .or. elementary_charge <= 0.0_real64) return
      energy_ev = 0.5_real64*mass*sum(velocity**2)
      if (.not. ieee_is_finite(energy_ev)) then
         energy_ev = 0.0_real64
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      energy_ev = energy_ev/elementary_charge
      if (.not. ieee_is_finite(energy_ev)) then
         energy_ev = 0.0_real64
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      ierr = SEEM_OK
   end subroutine compute_kinetic_energy_ev

   ! 将物理发射产额换算成宏粒子数的期望值：产额乘入射宏粒子权重，再除以产物权重。
   ! 该结果通常不是整数；后续可交由泊松采样或其他计数逻辑生成实际宏粒子数。
   subroutine compute_see_macro_mean(yield_value, incoming_weight, product_weight, mean, ierr)
      real(real64), intent(in) :: yield_value, incoming_weight, product_weight
      real(real64), intent(out) :: mean
      integer, intent(out) :: ierr

      mean = 0.0_real64
      ierr = SEEM_ERR_INVALID
      if (.not. ieee_is_finite(yield_value) .or. .not. ieee_is_finite(incoming_weight) .or. &
          .not. ieee_is_finite(product_weight)) return
      if (yield_value < 0.0_real64 .or. incoming_weight < 0.0_real64 .or. product_weight <= 0.0_real64) return

      mean = yield_value * incoming_weight
      if (.not. ieee_is_finite(mean)) then
         mean = 0.0_real64
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      mean = mean / product_weight
      if (.not. ieee_is_finite(mean)) then
         mean = 0.0_real64
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      ierr = SEEM_OK
   end subroutine compute_see_macro_mean

   ! 根据面电流密度估计当前步的发射宏粒子数，并用残差保存不足一个宏粒子的分数。
   ! expected = j*area*dt/(e*|Z|*macro_weight)；累加残差后取整数部分，余数留给下一步。
   subroutine accumulate_surface_emission(j, area, dt, elementary_charge, charge_number, &
         macro_weight, residual_in, n_emit, residual_out, ierr)
      real(real64), intent(in) :: j, area, dt, elementary_charge
      real(real64), intent(in) :: charge_number, macro_weight, residual_in
      integer, intent(out) :: n_emit, ierr
      real(real64), intent(out) :: residual_out
      real(real64) :: expected, total

      n_emit = 0
      residual_out = 0.0_real64
      ierr = SEEM_ERR_INVALID

      if (.not. ieee_is_finite(j) .or. .not. ieee_is_finite(area) .or. &
          .not. ieee_is_finite(dt) .or. .not. ieee_is_finite(elementary_charge) .or. &
          .not. ieee_is_finite(charge_number) .or. .not. ieee_is_finite(macro_weight) .or. &
          .not. ieee_is_finite(residual_in)) return
      if (j < 0.0_real64 .or. area <= 0.0_real64 .or. dt <= 0.0_real64 .or. &
          elementary_charge <= 0.0_real64 .or. abs(charge_number) <= 0.0_real64 .or. &
          macro_weight <= 0.0_real64 .or. residual_in < 0.0_real64 .or. &
          residual_in >= 1.0_real64) return

      ! 即使本步电流为零，也原样传递残差，避免跨时间步累计的分数被丢弃。
      residual_out = residual_in
      if (j <= 0.0_real64) then
         ierr = SEEM_OK
         return
      end if

      expected = j * area
      if (.not. ieee_is_finite(expected)) then
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      expected = expected * dt
      if (.not. ieee_is_finite(expected)) then
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      expected = expected / (elementary_charge * abs(charge_number))
      if (.not. ieee_is_finite(expected)) then
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      expected = expected / macro_weight
      if (.not. ieee_is_finite(expected)) then
         ierr = SEEM_ERR_OVERFLOW
         return
      end if

      total = expected + residual_in
      if (.not. ieee_is_finite(total) .or. total > real(huge(n_emit), real64)) then
         residual_out = 0.0_real64
         ierr = SEEM_ERR_OVERFLOW
         return
      end if

      ! floor 保证本步发射整数个宏粒子；total 的小数部分维持在 [0,1) 供下步使用。
      n_emit = floor(total)
      residual_out = total - real(n_emit, real64)
      ierr = SEEM_OK
   end subroutine accumulate_surface_emission

   ! 在线性递增的能量表中查找二分区间，并对相邻产额做线性插值。
   ! 表格无效时返回 SEEM_ERR_INVALID；能量越界时返回最近端点的产额值及 SEEM_OUT_OF_RANGE。
   subroutine interpolate_see_yield(energy_ev, table_energy, table_yield, yield_value, ierr)
      real(real64), intent(in) :: energy_ev
      real(real64), intent(in) :: table_energy(:), table_yield(:)
      real(real64), intent(out) :: yield_value
      integer, intent(out) :: ierr
      integer :: i, lower, upper, middle, n
      real(real64) :: fraction

      yield_value = 0.0_real64
      ierr = SEEM_ERR_INVALID
      n = size(table_energy)
      if (n < 2 .or. size(table_yield) /= n) return
      if (.not. ieee_is_finite(energy_ev)) return

      do i = 1, n
         if (.not. ieee_is_finite(table_energy(i)) .or. &
             .not. ieee_is_finite(table_yield(i))) return
         if (table_yield(i) < 0.0_real64) return
      end do
      do i = 2, n
         if (table_energy(i) <= table_energy(i-1)) return
      end do

      if (energy_ev < table_energy(1)) then
         yield_value = table_yield(1)
         ierr = SEEM_OUT_OF_RANGE
         return
      end if
      if (energy_ev > table_energy(n)) then
         yield_value = table_yield(n)
         ierr = SEEM_OUT_OF_RANGE
         return
      end if
      if (energy_ev <= table_energy(1)) then
         yield_value = table_yield(1)
         ierr = SEEM_OK
         return
      end if
      if (energy_ev >= table_energy(n)) then
         yield_value = table_yield(n)
         ierr = SEEM_OK
         return
      end if

      ! 二分搜索将区间逐次减半，避免对较长的产额表逐项线性扫描。
      lower = 1
      upper = n
      do while (upper - lower > 1)
         middle = lower + (upper - lower) / 2
         if (energy_ev < table_energy(middle)) then
            upper = middle
         else
            lower = middle
         end if
      end do

      fraction = (energy_ev - table_energy(lower)) / &
                 (table_energy(upper) - table_energy(lower))
      yield_value = table_yield(lower) + &
                    fraction * (table_yield(upper) - table_yield(lower))
      ierr = SEEM_OK
   end subroutine interpolate_see_yield

   ! 解析表面电流发射定义，字段依次为：边界组、物种、电流密度、温度和速度分布名称。
   ! 解析成功前先检查字段数、目标字符串容量及数值有效性，防止截断或无效参数进入后续初始化。
   subroutine parse_surface_current_definition(definition, group_name, species_name, &
         current_density, temperature, vdf_name, ierr)
      character(len=*), intent(in) :: definition
      character(len=*), intent(out) :: group_name, species_name, vdf_name
      real(real64), intent(out) :: current_density, temperature
      integer, intent(out) :: ierr
      character(len=512) :: tokens(6)
      integer :: token_count, ios

      group_name = ''
      species_name = ''
      vdf_name = ''
      current_density = 0.0_real64
      temperature = 0.0_real64
      ierr = SEEM_ERR_INVALID

      call tokenize_emission_definition(definition, tokens, token_count, ierr)
      if (ierr /= SEEM_OK) return
      ierr = SEEM_ERR_INVALID
      if (token_count /= 5) then
         ierr = SEEM_ERR_INVALID
         return
      end if
      if (len_trim(tokens(1)) > len(group_name) .or. &
          len_trim(tokens(2)) > len(species_name) .or. &
          len_trim(tokens(5)) > len(vdf_name)) then
         ierr = SEEM_ERR_INVALID
         return
      end if

      read(tokens(3), *, iostat=ios) current_density
      if (ios /= 0) return
      read(tokens(4), *, iostat=ios) temperature
      if (ios /= 0) return
      if (.not. ieee_is_finite(current_density) .or. .not. ieee_is_finite(temperature)) return
      if (current_density < 0.0_real64 .or. temperature <= 0.0_real64) return

      group_name = trim(tokens(1))
      species_name = trim(tokens(2))
      vdf_name = trim(tokens(5))
      if (len_trim(group_name) == 0 .or. len_trim(species_name) == 0 .or. len_trim(vdf_name) == 0) then
         group_name = ''
         species_name = ''
         vdf_name = ''
         return
      end if
      ierr = SEEM_OK
   end subroutine parse_surface_current_definition

   ! 解析二次电子发射定义，字段依次为：边界组、入射物种、电子物种、产额表文件和温度。
   ! 温度必须为有限正数；任一必需名称为空或超出输出缓冲区容量时，定义均视为无效。
   subroutine parse_secondary_emission_definition(definition, group_name, incident_species, &
         electron_species, table_file, temperature, ierr)
      character(len=*), intent(in) :: definition
      character(len=*), intent(out) :: group_name, incident_species, electron_species, table_file
      real(real64), intent(out) :: temperature
      integer, intent(out) :: ierr
      character(len=512) :: tokens(6)
      integer :: token_count, ios

      group_name = ''
      incident_species = ''
      electron_species = ''
      table_file = ''
      temperature = 0.0_real64
      ierr = SEEM_ERR_INVALID

      call tokenize_emission_definition(definition, tokens, token_count, ierr)
      if (ierr /= SEEM_OK) return
      ierr = SEEM_ERR_INVALID
      if (token_count /= 5) then
         ierr = SEEM_ERR_INVALID
         return
      end if
      if (len_trim(tokens(1)) > len(group_name) .or. &
          len_trim(tokens(2)) > len(incident_species) .or. &
          len_trim(tokens(3)) > len(electron_species) .or. &
          len_trim(tokens(4)) > len(table_file)) then
         ierr = SEEM_ERR_INVALID
         return
      end if

      read(tokens(5), *, iostat=ios) temperature
      if (ios /= 0) return
      if (.not. ieee_is_finite(temperature) .or. temperature <= 0.0_real64) return

      group_name = trim(tokens(1))
      incident_species = trim(tokens(2))
      electron_species = trim(tokens(3))
      table_file = trim(tokens(4))
      if (len_trim(group_name) == 0 .or. len_trim(incident_species) == 0 .or. &
          len_trim(electron_species) == 0 .or. len_trim(table_file) == 0) then
         group_name = ''
         incident_species = ''
         electron_species = ''
         table_file = ''
         return
      end if
      ierr = SEEM_OK
   end subroutine parse_secondary_emission_definition

   ! 将定义字符串拆成有界数量的字段：空格/制表符分隔，支持单引号或双引号包住字段。
   ! 感叹号开始行尾注释；字段过多、过长或引号未闭合时返回 SEEM_ERR_INVALID。
   subroutine tokenize_emission_definition(definition, tokens, token_count, ierr)
      character(len=*), intent(in) :: definition
      character(len=*), intent(out) :: tokens(:)
      integer, intent(out) :: token_count, ierr
      integer :: i, last, token_length
      character :: quote

      tokens = ''
      token_count = 0
      ierr = SEEM_ERR_INVALID
      last = len_trim(definition)
      i = 1

      do while (i <= last)
         do while (i <= last)
            if (definition(i:i) /= ' ' .and. definition(i:i) /= achar(9)) exit
            i = i + 1
         end do
         if (i > last) exit
         if (definition(i:i) == '!') exit

         token_count = token_count + 1
         if (token_count > size(tokens)) return
         token_length = 0
         quote = ' '
         if (definition(i:i) == '"' .or. definition(i:i) == "'") then
            quote = definition(i:i)
            i = i + 1
            do while (i <= last)
               if (definition(i:i) == quote) exit
               token_length = token_length + 1
               if (token_length > len(tokens(token_count))) return
               tokens(token_count)(token_length:token_length) = definition(i:i)
               i = i + 1
            end do
            if (i > last) return
            i = i + 1
            if (i <= last) then
               if (definition(i:i) /= ' ' .and. definition(i:i) /= achar(9) .and. &
                   definition(i:i) /= '!') return
            end if
         else
            do while (i <= last)
               if (definition(i:i) == ' ' .or. definition(i:i) == achar(9) .or. &
                   definition(i:i) == '!') exit
               token_length = token_length + 1
               if (token_length > len(tokens(token_count))) return
               tokens(token_count)(token_length:token_length) = definition(i:i)
               i = i + 1
            end do
         end if
      end do

      ierr = SEEM_OK
   end subroutine tokenize_emission_definition

   ! 两遍读取产额表：第一遍验证行格式并统计数据行数，第二遍分配数组并装入数值。
   ! 只接受至少两行、能量严格递增且产额非负的表；任何读取或校验失败都会关闭文件并释放已分配数组。
   subroutine read_see_yield_table(filename, table_energy, table_yield, ierr)
      character(len=*), intent(in) :: filename
      real(real64), allocatable, intent(out) :: table_energy(:), table_yield(:)
      integer, intent(out) :: ierr
      integer :: unit, ios, row_count, row, parse_status
      real(real64) :: energy, yield_value, previous_energy
      character(len=1024) :: line
      logical :: has_data

      ierr = SEEM_ERR_INVALID
      if (len_trim(filename) == 0) return
      open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) return

      ! 先统计有效数据行，才能一次性按精确尺寸分配输出数组。
      row_count = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios < 0) exit
         if (ios > 0) then
            close(unit)
            return
         end if
         call parse_see_yield_line(line, energy, yield_value, has_data, parse_status)
         if (parse_status /= SEEM_OK) then
            close(unit)
            return
         end if
         if (has_data) row_count = row_count + 1
      end do
      if (row_count < 2) then
         close(unit)
         return
      end if

      ! 回到文件开头后重读并填充数组；第二遍同时再次验证排序和行数一致性。
      rewind(unit, iostat=ios)
      if (ios /= 0) then
         close(unit)
         return
      end if
      allocate(table_energy(row_count), table_yield(row_count))
      row = 0
      previous_energy = -huge(previous_energy)
      do
         read(unit, '(A)', iostat=ios) line
         if (ios < 0) exit
         if (ios > 0) then
            close(unit)
            deallocate(table_energy, table_yield)
            return
         end if
         call parse_see_yield_line(line, energy, yield_value, has_data, parse_status)
         if (parse_status /= SEEM_OK) then
            close(unit)
            deallocate(table_energy, table_yield)
            return
         end if
         if (.not. has_data) cycle
         row = row + 1
         if (row > row_count) then
            close(unit)
            deallocate(table_energy, table_yield)
            return
         end if
         if (row > 1 .and. energy <= previous_energy) then
            close(unit)
            deallocate(table_energy, table_yield)
            return
         end if
         table_energy(row) = energy
         table_yield(row) = yield_value
         previous_energy = energy
      end do
      close(unit)
      if (row /= row_count) then
         deallocate(table_energy, table_yield)
         return
      end if

      ierr = SEEM_OK
   end subroutine read_see_yield_table

   ! 解析一行产额数据；空行和注释行不是错误，而是没有数据的有效行。
   ! 数据行必须恰含两个数值字段（能量和产额），并且二者有限、产额非负。
   subroutine parse_see_yield_line(raw_line, energy, yield_value, has_data, ierr)
      character(len=*), intent(in) :: raw_line
      real(real64), intent(out) :: energy, yield_value
      logical, intent(out) :: has_data
      integer, intent(out) :: ierr
      character(len=len(raw_line)) :: line
      character(len=1) :: extra
      integer :: comment_at, ios

      energy = 0.0_real64
      yield_value = 0.0_real64
      has_data = .false.
      ierr = SEEM_ERR_INVALID
      line = raw_line
      comment_at = index(line, '!')
      if (comment_at > 0) line(comment_at:) = ' '
      if (len_trim(adjustl(line)) == 0) then
         ierr = SEEM_OK
         return
      end if

      read(line, *, iostat=ios) energy, yield_value
      if (ios /= 0) return
      read(line, *, iostat=ios) energy, yield_value, extra
      if (ios >= 0) return
      if (.not. ieee_is_finite(energy) .or. .not. ieee_is_finite(yield_value)) return
      if (yield_value < 0.0_real64) return

      has_data = .true.
      ierr = SEEM_OK
   end subroutine parse_see_yield_line

   ! 从均值为 mean 的泊松分布中抽取非负整数计数，并在整数类型可能溢出时提前拒绝。
   ! 小均值使用概率递推反演；较大均值使用 PTRS 变换拒绝采样，提高高计数区间的效率。
   subroutine sample_poisson_count(mean, n, ierr)
      real(real64), intent(in) :: mean
      integer, intent(out) :: n, ierr
      real(real64) :: uniform, probability, cumulative
      real(real64) :: b, a, inverse_alpha, v_r, u, v, u_s, candidate, count_limit
      integer :: k

      n = 0
      ierr = SEEM_ERR_INVALID
      if (.not. ieee_is_finite(mean)) return
      if (mean < 0.0_real64) return
      count_limit = real(huge(n), real64)
      if (mean > count_limit) then
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      if (mean > count_limit - 10.0_real64 * sqrt(mean)) then
         ierr = SEEM_ERR_OVERFLOW
         return
      end if
      if (mean <= 0.0_real64) then
         ierr = SEEM_OK
         return
      end if

      if (mean < 10.0_real64) then
         ! 从 P(0)=exp(-mean) 开始递推 P(k)，累加到超过均匀随机数时即得到泊松样本。
         uniform = genrand64_real3()
         probability = exp(-mean)
         cumulative = probability
         k = 0
         do while (uniform > cumulative)
            k = k + 1
            probability = probability * mean / real(k, real64)
            cumulative = cumulative + probability
         end do
         n = k
         ierr = SEEM_OK
         return
      end if

      ! PTRS transformed rejection sampler (Hörmann) for moderate and large means.
      ! 快速接受条件减少常见情形的计算；其余候选值通过对数概率比较完成拒绝检验。
      b = 0.931_real64 + 2.53_real64 * sqrt(mean)
      a = -0.059_real64 + 0.02483_real64 * b
      inverse_alpha = 1.1239_real64 + 1.1328_real64 / (b - 3.4_real64)
      v_r = 0.9277_real64 - 3.6224_real64 / (b - 2.0_real64)

      do
         u = genrand64_real3() - 0.5_real64
         v = genrand64_real3()
         u_s = 0.5_real64 - abs(u)
         candidate = (2.0_real64 * a / u_s + b) * u + mean + 0.43_real64
         if (candidate < 0.0_real64) cycle
         if (candidate > count_limit) then
            ierr = SEEM_ERR_OVERFLOW
            return
         end if
         k = floor(candidate)

         if (u_s >= 0.07_real64 .and. v <= v_r) then
            n = k
            ierr = SEEM_OK
            return
         end if
         if (u_s < 0.013_real64 .and. v > u_s) cycle
         if (log(v * inverse_alpha / (a / (u_s * u_s) + b)) <= &
             -mean + real(k, real64) * log(mean) - log_gamma(real(k, real64) + 1.0_real64)) then
            n = k
            ierr = SEEM_OK
            return
         end if
      end do
   end subroutine sample_poisson_count

end module surface_emission_utils
