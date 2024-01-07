#! /usr/bin/env python3

import boto3
import argparse

class BigBadClass:
    def __init__(self, region_name='us-west-2'):
        self.rds_client = boto3.client('rds', region_name=region_name)

    def rds_create_modify_parameter_group(self, param_group_name, description, dbtype='postgres'):
        param_group_family='aurora-postgresql15'
        ssl_param_name = 'force_ssl'
        if dbtype == 'mysql':
            param_group_family="aurora-mysql8.0"
            ssl_param_name = 'require_secure_transport'
        try:
            # List all parameter groups
            response = self.rds_client.describe_db_parameter_groups()

            # Check if the specific parameter group exists
            exists = any(pg['DBParameterGroupName'] == param_group_name for pg in response['DBParameterGroups'])
            if not exists:
                # Create the parameter group
                self.rds_client.create_db_parameter_group(
                    DBParameterGroupName=param_group_name,
                    DBParameterGroupFamily=param_group_family,
                    Description=description
                )

            # Modify the parameter group to set 'require_secure_transport'
            self.rds_client.modify_db_parameter_group(
                DBParameterGroupName=param_group_name,
                Parameters=[
                    {
                        'ParameterName': ssl_param_name,
                        'ParameterValue': 'ON',
                        'ApplyMethod': 'immediate'
                    }
                ]
            
            )
        except Exception as e:
            print(f"Error in parameter group operation: {e}")


    def rds_create_serverless_v2_cluster(self, db_cluster_id, param_group, muser, mpass, dbtype='postgresql'):
        engineversion='15.1'
        if dbtype == 'mysql':
            engineversion='5.7.44'
     
        try:
            response = self.rds_client.create_db_cluster(
                DatabaseName='AWSEB',
                DBClusterIdentifier=db_cluster_id,
                Engine=f'aurora-{dbtype}',
                EngineMode='serverless',
                EngineVersion=engineversion,  # '8.0.35' or '5.7.mysql_aurora.2.10.0',
                MasterUsername=muser,
                MasterUserPassword=mpass,
                #DBParameterGroupName=param_group,
                ScalingConfiguration={
                    'AutoPause': True,
                    'MinCapacity': 1,
                    'MaxCapacity': 6
                },
                StorageEncrypted=True,
                EnablePerformanceInsights=False,
                EnableHttpEndpoint=True,
                DeletionProtection=False
            )
            return response
        except Exception as e:
            print(f"Error creating Aurora Serverless V2 MySQL cluster: {e}")


    def delete_aurora_db_cluster(self, db_cluster_identifier):
        try:
            response = self.rds_client.delete_db_cluster(
                DBClusterIdentifier=db_cluster_identifier,
                SkipFinalSnapshot=True  # Set to False if you want to take a final snapshot
            )
            return response
        except Exception as e:
            print(f"Error deleting Aurora DB cluster: {e}")

    def list_aurora_db_clusters(self):
        try:
            response = self.rds_client.describe_db_clusters()
            for cluster in response['DBClusters']:
                print(f"Cluster Identifier: {cluster['DBClusterIdentifier']}")
                print(f"Engine Version: {cluster['EngineVersion']}")
                print(f"Endpoint: {cluster.get('Endpoint', 'N/A')}")
                print(f"Port: {cluster['Port']}\n")
        except Exception as e:
            print(f"Error listing Aurora DB clusters: {e}")

def main():
    parser = argparse.ArgumentParser(description='Manage Aurora Serverless Clusters')
    subparsers = parser.add_subparsers(dest='command')

    # Create cluster parser
    create_parser = subparsers.add_parser('create', help='Create a new Aurora Serverless cluster')
    create_parser.add_argument('--dbtype', required=False, default='postgresql', help='DB type (postgres or mysql)')
    create_parser.add_argument('--identifier', required=True, help='DB Cluster Identifier')
    create_parser.add_argument('--username', required=True, help='Master Username')
    create_parser.add_argument('--password', required=True, help='Master Password')

    # Delete cluster parser
    delete_parser = subparsers.add_parser('delete', help='Delete an Aurora Serverless cluster')
    delete_parser.add_argument('--identifier', required=True, help='DB Cluster Identifier')

    # List clusters parser
    subparsers.add_parser('list', help='List all Aurora Serverless clusters')

    args = parser.parse_args()

    db_creator = BigBadClass()

    if args.command == 'create':
        #param_group_name = 'aurora-serverless-v2-param-group'        
        #param_group_description = 'Custom parameter group for Aurora Serverless V2 with enforced TLS'
        #db_creator.rds_create_modify_parameter_group(f'{param_group_name}-{args.dbtype}', param_group_description, dbtype=args.dbtype)
        param_group_name = None
        print(db_creator.rds_create_serverless_v2_cluster(args.identifier, param_group_name, args.username, args.password, args.dbtype))

    elif args.command == 'delete':
        print(db_creator.delete_aurora_db_cluster(args.identifier))

    elif args.command == 'list':
        db_creator.list_aurora_db_clusters()

if __name__ == "__main__":
    main()



