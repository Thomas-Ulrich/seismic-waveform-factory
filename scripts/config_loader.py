import os
import warnings
import yaml


class ConfigLoader:
    def __init__(self, config_file, schema):
        self.schema = schema
        self.config_file = config_file
        self.config = self._load_and_validate(config_file, schema)

    def _load_and_validate(self, config_file, schema):
        assert os.path.isfile(config_file), f"{config_file} not found"

        with open(config_file) as f:
            user_cfg = yaml.safe_load(f) or {}

        def check_type(section_name, value, expected_type):
            # fix: allow int when float expected
            if value is not None:
                if expected_type == float:
                    if not isinstance(value, (float, int)):
                        raise TypeError(
                            f"{section_name} expected float, got {type(value).__name__}"
                        )
                    value = float(value)  # coerce int → float
                # this is for t_before and t_after
                elif isinstance(expected_type, tuple):
                    raise_error = False
                    if expected_type == (float, dict):
                        if (not isinstance(value, expected_type)) and (
                            not isinstance(value, (int, dict))
                        ):
                            raise_error = True
                    else:
                        if not isinstance(value, expected_type):
                            raise_error = True
                    if raise_error:
                        raise TypeError(
                            f"{section_name} expected one of "
                            f"{[t.__name__ for t in expected_type]}, "
                            f"got {type(value).__name__}"
                        )
                elif not isinstance(value, expected_type):
                    raise TypeError(
                        f"{section_name} expected {expected_type.__name__}, "
                        f"got {type(value).__name__}"
                    )
            return value

        def when_condition_met(when_rules, context):
            """
            Evaluate 'when' condition dict (e.g. {'type': 'instaseis'})
            against the current context dictionary (the enclosing section).
            """
            for key, required_value in when_rules.items():
                current_value = context.get(key)
                if isinstance(required_value, (list, tuple)):
                    if current_value not in required_value:
                        return False
                else:
                    if current_value != required_value:
                        return False
            return True

        def validate_section(section_name, schema, data):
            expected_type = schema.get("type")
            sub_schema = schema.get("schema")
            # Fill defaults
            if data is None:
                if expected_type == dict:
                    data = {}
                elif expected_type == list:
                    data = []
                else:
                    data = schema.get("default")

            if expected_type == dict:
                if not isinstance(data, dict):
                    raise TypeError(
                        f"{section_name} expected dict, got {type(data).__name__}"
                    )
                validated = {}

                # --- warn for unknown keys ---
                unknown_keys = set(data.keys()) - set(sub_schema.keys())
                if unknown_keys:
                    warnings.warn(
                        f"Section {section_name} has unknown parameter(s): "
                        f"{', '.join(unknown_keys)}"
                    )

                for key, rules in sub_schema.items():
                    when_rules = rules.get("when")
                    if when_rules and not when_condition_met(when_rules, data):
                        # skip validation — not applicable in this context
                        continue

                    value = data.get(key, rules.get("default"))
                    validated[key] = validate_section(
                        f"{section_name}.{key}", rules, value
                    )
                return validated

            elif expected_type == list:
                if not isinstance(data, list):
                    raise TypeError(
                        f"{section_name} expected list, got {type(data).__name__}"
                    )
                validated = []
                if sub_schema:  # list of dicts
                    for i, item in enumerate(data):
                        # pass sub_schema as schema for each item
                        validated.append(
                            validate_section(
                                f"{section_name}[{i}]",
                                {"type": dict, "schema": sub_schema},
                                item,
                            )
                        )
                elif "element_type" in schema:  # list of scalars
                    elem_type = schema["element_type"]
                    for i, item in enumerate(data):
                        item = check_type(section_name[i], item, elem_type)
                        validated.append(item)
                else:
                    validated = data
                return validated

            else:
                # scalar
                value = data
                if value is None:
                    value = schema.get("default")

                value = check_type(section_name, value, expected_type)
                choices = schema.get("choices")
                if choices is not None and value not in choices:
                    raise ValueError(
                        f"{section_name} value {value} not in allowed choices {choices}"
                    )
                return value

        validated_cfg = {}
        for section, section_schema in schema.items():
            section_data = user_cfg.get(section, None)
            validated_cfg[section] = validate_section(
                section, section_schema, section_data
            )

        return validated_cfg

    def __getitem__(self, key):
        return self.config[key]

    def get(self, section, key, default=None):
        return self.config.get(section, {}).get(key, default)
