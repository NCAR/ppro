# UFO/ppro 代码分层分析

## 🎯 设计原则：参考 CRTM

### **分层标准：**

| 层次 | 职责 | 依赖 | 位置 |
|------|------|------|------|
| **核心物理层** | 纯算法，无框架依赖 | 无外部依赖 | ppro-lib |
| **UFO 接口层** | JEDI 框架集成 | JEDI/UFO 依赖 | ufo/ppro |

---

## 📊 **当前 ufo_PPRO_mod.F90 代码分析**

### ✅ **应该保留在 UFO 的代码（第 1-418 行）**

#### **1. ufo_PPRO_setup** (第 50-208 行)
```fortran
subroutine ufo_PPRO_setup(self, yaml_conf)
  use fckit_configuration_module  ! ← JEDI 依赖
  use oops_variables_mod          ! ← OOPS 依赖
  use obs_variables_mod           ! ← UFO 依赖
  
  ! YAML 配置解析
  call yaml_conf%get_or_die("microphysics option", micro_option)
  call yaml_conf%get_or_die("VertCoord", coord_name)
  
  ! GeoVaLs 变量设置
  call self%geovars%push_back(geovars_list)
  ...
end subroutine
```
**原因：** 纯 JEDI/UFO 框架代码，必须保留

---

#### **2. ufo_PPRO_simobs** (第 212-417 行)
```fortran
subroutine ufo_PPRO_simobs(self, geovals, obss, nvars, nlocs, hofx)
  use ufo_geovals_mod         ! ← UFO 依赖
  use obsspace_mod            ! ← IODA 依赖
  use vert_interp_mod         ! ← UFO 依赖
  
  ! 从 GeoVaLs 获取数据
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)
  
  ! 垂直插值
  call vert_interp_weights(...)
  call vert_interp_apply(...)
  
  ! ObsSpace 接口
  call obsspace_get_db(obss, "MetaData", "height", obsvcoord)
  
  ! 循环调用核心算子
  do iobs=1,nlocs
    call ufo_PPRO_sim1obs(...)
  enddo
end subroutine
```
**原因：** UFO 框架集成代码，必须保留

---

### ❓ **有争议的代码：ufo_PPRO_sim1obs** (第 420-770 行)

这个函数**既有 UFO 特定代码，也有纯物理逻辑**

#### **分析：**

```fortran
subroutine ufo_PPRO_sim1obs(iband, density_air, temp_air, qr, qs, qg, ...)
  use dualpol_op_tlad_mod ! for set trajectory
  use mpi                  ! ← 不应该在这里！
  
  ! 物理参数设置（第 480-489 行）
  qx_min = 0.0d0 
  ntx_min = 1.0D-3
  dmmax = 25.0d0
  ...
  
  ! 系数选择（第 493-509 行）
  if (iband == 1) then
    rain_coefs => sband_rain_coefs
    ...
  else
    rain_coefs => cband_rain_coefs
    ...
  endif
  
  ! 输入规整化（第 511-519 行）
  qrreg = max(qx_min(1), qr)
  ...
  
  ! 融化层方案选择（第 531-554 行）
  if (present(nr) .and. present(ns) .and. ...) then
    call melting_scheme_zhang24(...)  ! NSSL
  elseif (present(nr) .and. present(ns) .and. ...) then
    call melting_scheme_zhang24(...)  ! TCWA2
  elseif (present(nr) .and. ...) then
    call melting_scheme_zhang24(...)  ! Thompson
  else
    call melting_scheme_zhang24(...)  ! WSM6
  endif
  
  ! 水成物处理（第 578-711 行）
  ! 1. Pure rain
  call watercontent(...)
  call dm_z_2moment(...) or dm_z_wsm6(...)
  call dualpol_op_rain(...)
  
  ! 2-7. 其他水成物（雪、霰、雹等）
  ...
  
  ! 总体计算（第 715-740 行）
  call dualpol_op_total(...)
  
  ! 轨迹设置（第 744-768 行）
  if (present(set_trajectory) .and. set_trajectory) then
    call populate_traj(...)  ! ← TL/AD 专用
  endif
  
end subroutine
```

---

## 🎯 **我的建议：移动部分代码到 ppro-lib**

### **应该移到 ppro-lib 的部分：**

#### ✅ **移到 ppro-lib：核心算子函数**

创建一个新的高层接口：`ppro_compute_dualpol`

```fortran
! ppro-lib/src/dualpol_op_mod.f90 中添加

!> @brief 高层接口：计算双偏振变量（一个观测点）
!>
!> 这是 ppro-lib 的主要接口函数，封装了所有水成物处理逻辑
subroutine ppro_compute_dualpol(iband, density_air, temp_air, &
                                qr, qs, qg, qh, &
                                nr, ns, ng, nh, &
                                vg, vh, &
                                zh, zdr, kdp, phv, &
                                params, set_trajectory)
  ! 输入
  integer, intent(in) :: iband
  real(kind=8), intent(in) :: density_air, temp_air
  real(kind=8), intent(in) :: qr, qs, qg
  real(kind=8), intent(in), optional :: qh
  real(kind=8), intent(in), optional :: nr, ns, ng, nh
  real(kind=8), intent(in), optional :: vg, vh
  
  ! 输出
  real(kind=8), intent(out) :: zh, zdr, kdp, phv
  
  ! 参数
  type(ppro_params), intent(in), optional :: params
  logical, intent(in), optional :: set_trajectory
  
  ! ===== 移动这些逻辑到这里 =====
  ! 1. 参数设置（qx_min, ntx_min, dmmax）
  ! 2. 系数选择（S-band vs C-band）
  ! 3. 输入规整化
  ! 4. 融化层方案选择
  ! 5. 7 种水成物处理
  ! 6. 总体变量计算
  ! 7. 轨迹设置
  ! =================================
  
end subroutine ppro_compute_dualpol
```

**移动的代码：** ufo_PPRO_sim1obs 的第 480-768 行（约 290 行）

**保留在 UFO：** 只保留框架集成和调用
```fortran
! ufo_PPRO_mod.F90 (简化后)
subroutine ufo_PPRO_simobs(self, geovals, obss, nvars, nlocs, hofx)
  ! UFO 框架代码
  call ufo_geovals_get_var(...)
  call vert_interp_weights(...)
  
  ! 只初始化一次
  call init_coefs_once()
  
  do iobs=1,nlocs
    ! 简单的调用
    call ppro_compute_dualpol(iband(iobs), rho, t, qr, qs, qg, &
                              zh, zdr, kdp, phv, nr=nr)
    
    ! UFO 特定的后处理
    zh_dBZ = 10.0*log10(zh)
    zdr_dB = 10.0*log10(zdr)
    hofx(:,iobs) = [zh_dBZ, zdr_dB, kdp, phv]
  enddo
end subroutine
```

---

### **应该保留在 UFO 的部分：**

#### ✅ **保留在 UFO：框架集成**

| 代码 | 原因 | 依赖 |
|------|------|------|
| `ufo_PPRO_setup` | YAML 解析 | fckit_configuration |
| `ufo_PPRO_simobs` 的框架部分 | GeoVaLs、ObsSpace、插值 | UFO 框架 |
| ObsPPRO.cc/.h | C++ 接口 | JEDI C++ |
| ObsPPRO.interface.F90 | C++-Fortran 桥接 | JEDI |
| 变量名映射（opvars_list） | UFO 特定 | obs_variables_mod |

---

## 📊 **对比：CRTM 的分层**

### **CRTM-lib 中：**
```fortran
! 核心 API
CRTM_Forward(atm, sfc, geo, chinfo, rts)
  ↓
完全独立的物理算法
```

### **UFO 中（ufo_radiancecrtm_mod.F90）：**
```fortran
subroutine ufo_radiancecrtm_simobs(...)
  ! UFO 框架代码
  call Load_Atm_Data(geovals, atm)     ! ← UFO 特定
  call Load_Sfc_Data(geovals, sfc)     ! ← UFO 特定
  call Load_Geom_Data(obss, geo)       ! ← UFO 特定
  
  ! 调用 CRTM 核心
  call CRTM_Forward(atm, sfc, geo, ...)  ! ← 纯物理
  
  ! UFO 后处理
  hofx(:,:) = rts%Brightness_Temperature  ! ← UFO 特定
end subroutine
```

---

## 💡 **我的建议：移动 ufo_PPRO_sim1obs 到 ppro-lib**

### **新的分层设计：**

```
┌─────────────────────────────────────────────────────────┐
│                    UFO 层                                │
├─────────────────────────────────────────────────────────┤
│ ufo_PPRO_simobs (ufo_PPRO_mod.F90)                     │
│  ├─ YAML 配置解析 ✅                                    │
│  ├─ GeoVaLs 数据提取 ✅                                 │
│  ├─ 垂直插值 ✅                                         │
│  ├─ ObsSpace 接口 ✅                                    │
│  └─ 循环调用核心算子                                    │
│       ↓                                                  │
├─────────────────────────────────────────────────────────┤
│              接口：ppro_compute_dualpol                 │
├─────────────────────────────────────────────────────────┤
│                  ppro-lib 层                            │
├─────────────────────────────────────────────────────────┤
│ ppro_compute_dualpol (dualpol_op_mod.f90)              │
│  ├─ 参数设置（qx_min, dmmax） ⬇️ 移动               │
│  ├─ 系数选择（S/C-band） ⬇️ 移动                     │
│  ├─ 融化层方案选择 ⬇️ 移动                           │
│  ├─ 7 种水成物处理 ⬇️ 移动                           │
│  ├─ dualpol_op_rain ✅ 已在                           │
│  ├─ dualpol_op_icephase ✅ 已在                       │
│  ├─ dualpol_op_total ✅ 已在                          │
│  └─ populate_traj ✅ 已在 TL/AD 模块                  │
└─────────────────────────────────────────────────────────┘
```

---

## 📋 **具体建议**

### **方案 A：激进移动（推荐，更像 CRTM）**

#### 移动到 ppro-lib：

**1. ufo_PPRO_sim1obs 的核心部分** (约 350 行)
- ✅ 物理参数设置
- ✅ 融化层方案选择逻辑
- ✅ 7 种水成物的完整处理流程
- ✅ 系数选择逻辑

**保留在 UFO：**
- ✅ YAML 配置解析（setup）
- ✅ GeoVaLs 和 ObsSpace 接口（simobs）
- ✅ 垂直插值
- ✅ 单位转换（dBZ, dB）
- ✅ 缺失值处理

#### **代码示例：**

**ppro-lib（新增）：**
```fortran
! ppro-lib/src/dualpol_op_mod.f90

subroutine ppro_compute_point(iband, scheme_type, &
                               density_air, temp_air, &
                               qr, qs, qg, qh, &
                               nr, ns, ng, nh, vg, vh, &
                               zh, zdr, kdp, phv, &
                               coef_path, set_trajectory)
  ! 所有输入
  integer, intent(in) :: iband  ! 1=S-band, 2=C-band
  character(len=*), intent(in) :: scheme_type  ! 'WSM6', 'Thompson', 'NSSL'
  real(kind=8), intent(in) :: density_air, temp_air
  real(kind=8), intent(in) :: qr, qs, qg
  real(kind=8), intent(in), optional :: qh, nr, ns, ng, nh, vg, vh
  
  ! 所有输出
  real(kind=8), intent(out) :: zh, zdr, kdp, phv
  
  ! 可选参数
  character(len=*), intent(in), optional :: coef_path
  logical, intent(in), optional :: set_trajectory
  
  ! === 这里包含所有 ufo_PPRO_sim1obs 的核心逻辑 ===
  ! 1. 初始化系数（只一次）
  ! 2. 选择 S/C-band 系数
  ! 3. 融化层方案
  ! 4. 7 种水成物处理
  ! 5. 总体计算
  ! 6. 轨迹设置
  
end subroutine ppro_compute_point
```

**UFO（简化后）：**
```fortran
! ufo_PPRO_mod.F90

subroutine ufo_PPRO_simobs(self, geovals, obss, nvars, nlocs, hofx)
  ! ... UFO 框架代码 ...
  
  ! 初始化系数（只一次）
  call ppro_init_coefs_once()
  
  do iobs=1,nlocs
    ! 从 fields 提取数据
    rho = fields(1,iobs)
    t = fields(2,iobs)
    qr = 1000.0*fields(5,iobs)
    ...
    
    ! 调用 ppro-lib 核心函数
    call ppro_compute_point(iband(iobs), trim(self%micro_option), &
                            rho, t, qr, qs, qg, qh, &
                            nr, ns, ng, nh, vg, vh, &
                            zh, zdr, kdp, phv, &
                            set_trajectory=.false.)
    
    ! UFO 后处理：转换单位
    if (zh > 0.0) then
      zh_dBZ = 10.0*log10(zh)
    else
      zh_dBZ = missing
    endif
    
    hofx(obsvarindex,iobs) = zh_dBZ
    ...
  enddo
end subroutine
```

**收益：**
- ✅ ppro-lib 完全独立，可单独测试
- ✅ UFO 代码更简洁（从 775 行减少到 ~400 行）
- ✅ 物理算法可以在其他项目中复用

---

### **方案 B：保守方案（当前状态）**

#### 保持现状：
- ppro-lib 只包含低层函数
- ufo_PPRO_sim1obs 保留在 UFO

**优点：**
- 改动最小
- 当前已经能工作

**缺点：**
- ppro-lib 不够独立
- 用户必须理解 UFO 代码才能使用 ppro-lib

---

## 🎯 **参考 CRTM 的完整度**

### CRTM-lib 提供的 API：

```fortran
! 完整的高层接口
CRTM_Forward(Atmosphere, Surface, Geometry, RTSolution)
  ↓
用户只需准备输入结构，CRTM 处理所有细节
```

### 如果 ppro-lib 也这样：

```fortran
! ppro-lib 的高层接口
ppro_compute_point(iband, scheme, density_air, temp, qr, qs, qg, ...)
  ↓
用户只需提供输入，ppro-lib 处理所有物理计算
```

---

## 📊 **代码行数对比**

| 层次 | 当前（方案 B） | 建议（方案 A） | 变化 |
|------|--------------|--------------|------|
| **ppro-lib** | 968 行 | ~1300 行 | +350 行 |
| **ufo/ppro (Fortran)** | 775 行 | ~400 行 | -375 行 |
| **总代码** | 1743 行 | 1700 行 | -43 行（更简洁） |

---

## 🎯 **我的最终建议**

### **推荐：方案 A（激进移动）**

**原因：**
1. ✅ **更像 CRTM** - ppro-lib 成为完整的独立库
2. ✅ **更易复用** - 其他项目可以直接使用 ppro-lib
3. ✅ **更易测试** - 可以独立测试核心算法
4. ✅ **职责分离** - UFO 只负责框架集成，ppro-lib 负责物理

**需要移动的函数/代码：**
```
ufo_PPRO_mod.F90 第 420-770 行
  → ppro-lib/src/dualpol_op_mod.f90
  → 重命名为 ppro_compute_point
```

---

## 🤔 **具体问题：**

### **明确应该移动的：**
1. ✅ 融化层方案选择逻辑（第 531-554 行）
2. ✅ 7 种水成物处理流程（第 578-711 行）
3. ✅ 物理参数设置（qx_min, dmmax）
4. ✅ S/C-band 系数选择

### **应该保留在 UFO 的：**
1. ✅ YAML 配置解析
2. ✅ GeoVaLs 插值
3. ✅ ObsSpace 接口
4. ✅ 单位转换（linear → dBZ）
5. ✅ 缺失值处理

---

## ❓ **您的选择？**

1. **方案 A**：移动 ufo_PPRO_sim1obs 到 ppro-lib（更彻底的模块化）
2. **方案 B**：保持当前状态（最小改动）
3. **方案 C**：只移动部分（折中）

我可以根据您的选择实施相应的重构！🚀

