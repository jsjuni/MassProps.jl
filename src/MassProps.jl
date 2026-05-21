"Calculate mass properties and their uncertainties for composite structures."
module MassProps

    using RollupTree
    using LinearAlgebra

    export get_mass_props, get_mass_props_unc, get_mass_props_and_unc,
            set_mass_props, set_mass_props_unc, set_mass_props_and_unc,
            combine_mass_props, combine_mass_props_unc, combine_mass_props_and_unc,
            set_poi_conv_plus, set_poi_conv_minus, set_poi_conv_from_target,
            update_mass_props, update_mass_props_unc, update_mass_props_and_unc,
            validate_mass_props, validate_mass_props_unc, validate_mass_props_and_unc,
            validate_mass_props_table, validate_mass_props_and_unc_table,
            rollup_mass_props, rollup_mass_props_unc, rollup_mass_props_and_unc,
            X, Y, Z

    const X = 1
    const Y = 2
    const Z = 3

    """
        get_mass_props(table, id)

        Returns a named tuple with the following fields:
    - `mass`: the mass of the component
    - `center_mass`: a 3-element vector containing the x, y, z coordinates of the center of mass
    - `inertia`: a 3x3 matrix containing the inertia tensor
    - `point`: a boolean indicating whether the object is a point mass
    """
    function get_mass_props(table, id)
        row = df_get_row_by_id(table, id)
        poi_factor = row.POIconv == "-" ? 1.0 : -1.0

        (
            mass = row.mass,

            center_mass = [row.Cx; row.Cy; row.Cz],

            inertia = [             row.Ixx poi_factor * row.Ixy poi_factor * row.Ixz;
                       poi_factor * row.Ixy              row.Iyy poi_factor * row.Iyz;
                       poi_factor * row.Ixz poi_factor * row.Iyz              row.Izz
            ],

           point = row.Ipoint
        )
    end

    """
        get_mass_props_unc(table, id)

        Returns a named tuple with the following fields:
    - `sigma_mass`: the uncertainty in the mass
    - `sigma_center_mass`: a 3-element vector containing the uncertainties in the x, y, z coordinates of the center of mass
    - `sigma_inertia`: a 3x3 matrix containing the uncertainties in the inertia tensor
    """
    function get_mass_props_unc(table, id)
        row = df_get_row_by_id(table, id)

        (
            sigma_mass = row.sigma_mass,

            sigma_center_mass = [row.sigma_Cx; row.sigma_Cy; row.sigma_Cz],

            sigma_inertia = [row.sigma_Ixx row.sigma_Ixy row.sigma_Ixz;
                             row.sigma_Ixy row.sigma_Iyy row.sigma_Iyz;
                             row.sigma_Ixz row.sigma_Iyz row.sigma_Izz
            ]
        )
    end

    get_mass_props_and_unc(table, id) = merge(get_mass_props(table, id), get_mass_props_unc(table, id))

    
    """
        set_mass_props(table, id, mp)

        Sets the mass properties for the given `id` in the `table` using the values from the named tuple `mp`, which should have the following fields:
    - `mass`: the mass of the component
    - `center_mass`: a 3-element vector containing the x, y, z coordinates of the center of mass
    - `inertia`: a 3x3 matrix containing the inertia tensor
    - `point`: a boolean indicating whether the object is a point mass
    - `poi_conv`: a string indicating the POI convention ("+" or "-")
    """
    function set_mass_props(table, id, mp)
        
        cm = mp.center_mass
        it = (mp.inertia + mp.inertia') / 2 # ensure symmetry

        poi_factor = mp.poi_conv == "+" ? -1.0 : (mp.poi_conv == "-" ? 1.0 : error("invalid POI convention"))

        values = (
            mass = mp.mass,

            Cx = cm[X],
            Cy = cm[Y],
            Cz = cm[Z],

            Ixx = it[X, X],
            Iyy = it[Y, Y],
            Izz = it[Z, Z],

            Ixy = poi_factor * it[X, Y],
            Ixz = poi_factor * it[X, Z],
            Iyz = poi_factor * it[Y, Z],

            Ipoint = mp.point,
            POIconv = mp.poi_conv
        )

        df_set_row_by_id(table, id, values)
    end

    """
        set_mass_props_unc(table, id, mp_unc)

        Sets the mass properties uncertainties for the given `id` in the `table` using the values from the named tuple `mp_unc`, which should have the following fields:
    - `sigma_mass`: the uncertainty in the mass
    - `sigma_center_mass`: a 3-element vector containing the uncertainties in the x, y, z coordinates of the center of mass
    - `sigma_inertia`: a 3x3 matrix containing the uncertainties in the inertia tensor
    """
    function set_mass_props_unc(table, id, mp_unc)

        sigma_it = (mp_unc.sigma_inertia + mp_unc.sigma_inertia') / 2 # ensure symmetry
        
        values = (
            sigma_mass = mp_unc.sigma_mass,

            sigma_Cx = mp_unc.sigma_center_mass[X],
            sigma_Cy = mp_unc.sigma_center_mass[Y],
            sigma_Cz = mp_unc.sigma_center_mass[Z],

            sigma_Ixx = sigma_it[X, X],
            sigma_Iyy = sigma_it[Y, Y],
            sigma_Izz = sigma_it[Z, Z],

            sigma_Ixy = sigma_it[X, Y],
            sigma_Ixz = sigma_it[X, Z],
            sigma_Iyz = sigma_it[Y, Z]
        )

        df_set_row_by_id(table, id, values)
    end

    """
        set_mass_props_and_unc(table, id, mpu)

        Sets both the mass properties and their uncertainties for the given `id` in the `table` using the values from the named tuple `mpu`, which should have the following fields:
    - `mass`: the mass of the component
    - `center_mass`: a 3-element vector containing the x, y, z coordinates of the center of mass
    - `inertia`: a 3x3 matrix containing the inertia tensor
    - `point`: a boolean indicating whether the object is a point mass
    - `sigma_mass`: the uncertainty in the mass
    - `sigma_center_mass`: a 3-element vector containing the uncertainties in the x, y, z coordinates of the center of mass
    - `sigma_inertia`: a 3x3 matrix containing the uncertainties in the inertia tensor
    - `poi_conv`: a string indicating the POI convention ("+" or "-")
    """
    set_mass_props_and_unc(table, id, mpu) = set_mass_props_unc(set_mass_props(table, id, mpu), id, mpu)

    """
        combine_mass_props(mpl)

    Combines a list of mass properties named tuples `mpl` into a single mass properties named tuple representing the combined properties of the components. The input list should contain named tuples with the following fields:
    - `mass`: the mass of the component
    - `center_mass`: a 3-element vector containing the x, y, z coordinates of the center of mass
    - `inertia`: a 3x3 matrix containing the inertia tensor
    - `point`: a boolean indicating whether the object is a point mass
    """
    function combine_mass_props(mpl)
        
        (
            mass = (mass = sum(map(mp -> mp.mass, mpl))),

            center_mass = (center_mass = sum(map(mp -> mp.mass .* mp.center_mass, mpl)) ./ mass),

            inertia = sum(
                map(mp -> begin
                    d = mp.center_mass - center_mass
                    Q = d .* d'
                    M = mp.mass .* (tr(Q) * I - Q)
                    mp.point ? M : mp.inertia .+ M
                end, mpl)
            ),

            point = false

        )

    end

    """
        combine_mass_props_unc(mpul, amp)

    Combines a list of mass properties uncertainties named tuples `mpul` into a single mass properties uncertainty named tuple representing the combined uncertainties of the components, given the combined mass properties `amp`. The input list should contain named tuples with the following fields:
    - `sigma_mass`: the uncertainty in the mass
    - `sigma_center_mass`: a 3-element vector containing the uncertainties in the x, y, z coordinates of the center of mass
    - `sigma_inertia`: a 3x3 matrix containing the uncertainties in the inertia tensor
    - `amp`: a named tuple containing the combined mass properties with fields `mass`, `center_mass`, and `inertia` as returned by `combine_mass_props`.
    """
    function combine_mass_props_unc(mpul, amp)

        return merge(amp,
            (
                sigma_mass = (sigma_mass = sqrt(sum(map(mpu -> mpu.sigma_mass^2, mpul)))),

                sigma_center_mass = (sigma_center_mass = sqrt.(
                    sum(map(mpu -> 
                        (mpu.mass .* mpu.sigma_center_mass).^2 .+
                        (mpu.sigma_mass .* (mpu.center_mass - amp.center_mass)).^2,
                        mpul
                    ))
                ) ./ amp.mass),

                sigma_inertia = sqrt.(
                    sum(map(mpu -> begin
                        d = mpu.center_mass - amp.center_mass
                        P = d .* mpu.sigma_center_mass'
                        p = diag(P)
                        Q = d .* d'

                        M1 = P  - diagm(p - 2 .* view(p, [Y, X, X]))
                        M2 = P' - diagm(p - 2 .* view(p, [Z, Z, Y]))
                        M3 = Q  - tr(Q) * I
                        M4 = mpu.mass^2 .* (M1.^2 .+ M2.^2) .+ (mpu.sigma_mass .* M3).^2
                        mpu.point ? M4 : mpu.sigma_inertia.^2 .+ M4
                    end,
                    mpul))
                )
            )
        )

    end

    """
        combine_mass_props_and_unc(mpul)

    Combines a list of mass properties named tuples with their uncertainties `mpul` into a single named tuple representing both the combined mass properties and their uncertainties. The input list should contain named tuples with the following fields:
    - `mass`: the mass of the component
    - `center_mass`: a 3-element vector containing the x, y, z coordinates of the center of mass
    - `inertia`: a 3x3 matrix containing the inertia tensor
    - `point`: a boolean indicating whether the object is a point mass
    - `sigma_mass`: the uncertainty in the mass
    - `sigma_center_mass`: a 3-element vector containing the uncertainties in the x, y, z coordinates of the center of mass
    - `sigma_inertia`: a 3x3 matrix containing the uncertainties in the inertia tensor
    """
    combine_mass_props_and_unc(mpul) = combine_mass_props_unc(mpul, combine_mass_props(mpul))

    """
        set_poi_conv_plus(df, target, mp)

        Sets the POI convention to "+" for the given `target` in the `df` when updating mass properties.
    """
    set_poi_conv_plus(df, target, mp) = merge(mp, (poi_conv = "+",))
    """
        set_poi_conv_minus(df, target, mp)

        Sets the POI convention to "-" for the given `target` in the `df` when updating mass properties.
    """
    set_poi_conv_minus(df, target, mp) = merge(mp, (poi_conv = "-",))

    """
        set_poi_conv_from_target(df, target, mp)

        Sets the POI convention based on the value in the `df` for the given `target` when updating mass properties. The function retrieves the POI convention from the `df` using `df_get_by_id` and merges it into the mass properties named tuple `mp`.
    """
    set_poi_conv_from_target(df, target, mp) = merge(mp, (poi_conv = df_get_by_id(df, target, :POIconv),))

    """
        update_mass_props(df, target, sources, override = set_poi_conv_from_target)

        Updates the mass properties for the `target` in the `df` by combining the mass properties of the `sources`. The `override` function is used to modify the combined mass properties before setting them in the `df`. By default, it uses `set_poi_conv_from_target` to determine the POI convention based on the target's existing value in the `df`.
    """
    function update_mass_props(df, target, sources, override = set_poi_conv_from_target)
        update_prop(
            df,
            target,
            sources,
            set_mass_props,
            get_mass_props,
            combine = combine_mass_props,
            override = override
        )
    end

    """
        update_mass_props_unc(df, target, sources)

        Updates the mass properties uncertainties for the `target` in the `df` by combining the mass properties uncertainties of the `sources`. The function uses `combine_mass_props_unc` to combine the uncertainties based on the combined mass properties of the sources and sets the result in the `df`.
    """
    function update_mass_props_unc(df, target, sources)
        update_prop(
            df,
            target,
            sources,
            set_mass_props_unc,
            get_mass_props_and_unc,
            combine = l -> combine_mass_props_unc(l, get_mass_props(df, target))
        )
    end

    """
        update_mass_props_and_unc(df, target, sources, override = set_poi_conv_from_target)

        Updates both the mass properties and their uncertainties for the `target` in the `df` by combining the mass properties and uncertainties of the `sources`. The `override` function is used to modify the combined mass properties and uncertainties before setting them in the `df`. By default, it uses `set_poi_conv_from_target` to determine the POI convention based on the target's existing value in the `df`.
    """
    function update_mass_props_and_unc(df, target, sources, override = set_poi_conv_from_target)
           update_prop(
            df,
            target,
            sources,
            set_mass_props_and_unc,
            get_mass_props_and_unc,
            combine = combine_mass_props_and_unc,
            override = override
        )
    end

    """
    validate_mass_props(mp)

        Validates the mass properties named tuple `mp` to ensure that it contains valid values for mass, center of mass, inertia tensor, and point mass flag. The function checks for missing or invalid values and throws errors if any issues are found. The expected fields in `mp` are:
    - `mass`: the mass of the component (must be a positive real number)
    - `center_mass`: a 3-element vector containing the x, y, z coordinates of the center of mass (must be a vector of real numbers with three components)
    - `inertia`: a 3x3 matrix containing the inertia tensor (must be a 3x3 matrix of real numbers that is positive definite and satisfies triangle inequalities)
    - `point`: a boolean indicating whether the object is a point mass (must be a boolean)
    """
    function validate_mass_props(mp)
        
        ismissing(mp.mass) && error("mass is missing")
        isnothing(mp.mass) && error("mass is nothing")
        mp.mass isa Real || error("mass must be a real number")
        mp.mass > 0.0 || error("mass must be positive")

        any(ismissing.(mp.center_mass)) && error("center mass component is missing")
        any(isnothing.(mp.center_mass)) && error("center mass component is nothing")
        mp.center_mass isa AbstractVector{<:Real} || error("center mass must be a vector of real numbers")
        length(mp.center_mass) == 3 || error("center mass must have three components")

        any(ismissing.(mp.inertia)) && error("inertia component is missing")
        any(isnothing.(mp.inertia)) && error("inertia component is nothing")
        mp.inertia isa AbstractMatrix{<:Real} || error("inertia must be a matrix of real numbers")
        size(mp.inertia) == (3, 3) || error("inertia must be a 3x3 matrix")

        ev = eigen(mp.inertia)
        all(isreal.(ev.values)) || error("inertia matrix must have real eigenvalues")
        any(ev.values .<= 0.0) && error("inertia matrix must be positive definite")
        all([
            ev.values[1] <= ev.values[2] + ev.values[3],
            ev.values[2] <= ev.values[1] + ev.values[3],
            ev.values[3] <= ev.values[1] + ev.values[2]
        ]) || error("inertia matrix must satisfy triangle inequalities")

        mp.point isa Bool || error("point must be a boolean")

        return true

    end

    """
        validate_mass_props_unc(mpu)

    Validates the mass properties uncertainty named tuple `mpu` to ensure that it contains valid values for mass uncertainty, center of mass uncertainty, inertia uncertainty, and point mass flag uncertainty. The function checks for missing or invalid values and throws errors if any issues are found. The expected fields in `mpu` are:
    - `sigma_mass`: the uncertainty in the mass of the component (must be a non-negative real number)
    - `sigma_center_mass`: a 3-element vector containing the x, y, z coordinates of the center of mass uncertainty (must be a vector of real numbers with three components)
    - `sigma_inertia`: a 3x3 matrix containing the inertia tensor uncertainty (must be a 3x3 matrix of real numbers that is positive definite and satisfies triangle inequalities)
    - `sigma_point`: a boolean indicating whether the object is a point mass (must be a boolean)
    """
    function validate_mass_props_unc(mpu)

        mpu.sigma_mass isa Real || error("mass uncertainty must be a real number")
        mpu.sigma_mass >= 0.0 || error("mass uncertainty must be non-negative")

        any(ismissing.(mpu.sigma_center_mass)) && error("center mass uncertainty component is missing")
        any(isnothing.(mpu.sigma_center_mass)) && error("center mass uncertainty component is nothing")
        mpu.sigma_center_mass isa AbstractVector{<:Real} || error("center mass uncertainty must be a vector of real numbers")
        length(mpu.sigma_center_mass) == 3 || error("center mass uncertainty must have three components")
        any(mpu.sigma_center_mass .< 0.0) && error("center mass uncertainty must be non-negative")

        any(ismissing.(mpu.sigma_inertia)) && error("inertia uncertainty component is missing")
        any(isnothing.(mpu.sigma_inertia)) && error("inertia uncertainty component is nothing")
        mpu.sigma_inertia isa AbstractMatrix{<:Real} || error("inertia uncertainty must be a matrix of real numbers")
        size(mpu.sigma_inertia) == (3, 3) || error("inertia uncertainty must be a 3x3 matrix")
        any(mpu.sigma_inertia .< 0.0) && error("inertia uncertainty must be non-negative")

        return true

    end

    """
        validate_mass_props_and_unc(mpu)

    Validates both the mass properties and their uncertainties in the named tuple `mpu` by calling `validate_mass_props` and `validate_mass_props_unc`. The function returns true if both validations pass, otherwise it throws an error with details about the issues found in either the mass properties or their uncertainties.
    """
    validate_mass_props_and_unc(mpu) = validate_mass_props(mpu) && validate_mass_props_unc(mpu)

    """
        validate_mass_props_table(tree, df)

    Validates the mass properties for all entries in the `df` using the `tree` structure to determine the order of validation. The function uses `validate_ds` to apply the `validate_mass_props` function to each entry in the `df` based on the IDs obtained from `df_get_ids` and the mass properties obtained from `get_mass_props`. If any entry fails validation, an error is thrown with details about the issues found.
    """
    validate_mass_props_table(tree, df) = validate_ds(tree, df, df_get_ids, get_mass_props, validate_mass_props)

    """
        validate_mass_props_and_unc_table(tree, df)

    Validates both the mass properties and their uncertainties for all entries in the `df` using the `tree` structure to determine the order of validation. The function uses `validate_ds` to apply the `validate_mass_props_and_unc` function to each entry in the `df` based on the IDs obtained from `df_get_ids` and the combined mass properties and uncertainties obtained from `get_mass_props_and_unc`. If any entry fails validation, an error is thrown with details about the issues found.
    """
    validate_mass_props_and_unc_table(tree, df) = validate_ds(tree, df, df_get_ids, get_mass_props_and_unc, validate_mass_props_and_unc)

    """
        rollup_mass_props(tree, df, validate_df = validate_mass_props_table)

    Rolls up the mass properties for all entries in the `df` using the `tree` structure to determine the order of rollup. The function uses `rollup` to apply the `update_mass_props` function to each entry in the `df` based on the IDs obtained from `df_get_ids`, the mass properties obtained from `get_mass_props`, and the validation function specified by `validate_df`. By default, it uses `validate_mass_props_table` to validate the mass properties after each update. If any entry fails validation during the rollup process, an error is thrown with details about the issues found.
    """
    rollup_mass_props(tree, df, validate_df = validate_mass_props_table) = rollup(tree, df, update_mass_props, validate_df)

    """
        rollup_mass_props_unc(tree, df, validate_df = validate_mass_props_and_unc_table)

    Rolls up the mass properties uncertainties for all entries in the `df` using the `tree` structure to determine the order of rollup. The function uses `rollup` to apply the `update_mass_props_unc` function to each entry in the `df` based on the IDs obtained from `df_get_ids`, the combined mass properties and uncertainties obtained from `get_mass_props_and_unc`, and the validation function specified by `validate_df`. By default, it uses `validate_mass_props_and_unc_table` to validate both the mass properties and their uncertainties after each update. If any entry fails validation during the rollup process, an error is thrown with details about the issues found.
    """
    rollup_mass_props_unc(tree, df, validate_df = validate_mass_props_and_unc_table) = rollup(tree, df, update_mass_props_unc, validate_df)

    """
        rollup_mass_props_and_unc(tree, df, validate_df = validate_mass_props_and_unc_table)

    Rolls up both the mass properties and their uncertainties for all entries in the `df` using the `tree` structure to determine the order of rollup. The function uses `rollup` to apply the `update_mass_props_and_unc` function to each entry in the `df` based on the IDs obtained from `df_get_ids`, the combined mass properties and uncertainties obtained from `get_mass_props_and_unc`, and the validation function specified by `validate_df`. By default, it uses `validate_mass_props_and_unc_table` to validate both the mass properties and their uncertainties after each update. If any entry fails validation during the rollup process, an error is thrown with details about the issues found.
    """
    rollup_mass_props_and_unc(tree, df, validate_df = validate_mass_props_and_unc_table) = rollup(tree, df, update_mass_props_and_unc, validate_df)

end